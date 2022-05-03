# conciliation tools
library(tidyverse)
library(DemoTools)

t      <- 2015:2020
pop    <- tibble(age=0:100,pop=sample(10e5:10e6,101),t=2015)
births <- tibble(yb=2015:2020, B=sample(10000:8000,6)) %>% 
          filter(yb %in% t)
deaths <- crossing(yd=1900:2020,age=0:100) %>% 
          group_by(yd,age) %>% 
          mutate(D=sample(1:100,1)) %>% 
          filter(yd %in% t)
nmigr <- crossing(ym=1900:2020,age=0:100) %>% 
          group_by(ym,age) %>% 
          mutate(NM = 0) %>% 
          filter(ym %in% t)
deaths <- read.csv("data/D.csv") %>% group_by(anio, EDAD) %>% summarise(D = sum(n)) %>%
          rename(yd = anio, age=EDAD)
births <- read.csv("data/B.csv") %>% group_by(anio) %>% summarise(B = sum(n)) %>% rename(yb=anio)
pop <- read.csv("data/pop.csv") %>% 
        filter(codgeo==1, date==2010.5, tipo=="base", sex=="T") %>% 
        select(age, pop=valor) %>% 
        mutate(t=2010)
nmigr = NULL
debugonce(lexis_pop)
pop_out <- lexis_pop(pop, births, deaths)[[1]]
pop_out %>% ggplot(aes(age, pop)) + geom_line() + facet_grid(~t)


# going further pop: all with same time interval
lexis_pop <- function(pop, births, deaths, nmigr = NULL){
  
  delta <- unique(diff(births$yb))
  
  # split death triangles: return yb, yd, age, D    
  deaths <- split_triangles_D(deaths, births)
  
  # split smigr triangles: return yb, yd, age, M    
  if(is.null(nmigr)){
    nmigr <- deaths %>% rename(ym=yd, NM=D) %>% mutate(NM=0)
  }else{
    nmigr <- split_triangles_NM(nmigr)  
  }

  # rename pop
  pop <- pop  %>% mutate(yb = t - age - delta)
  
  # new cohorts 
  pop_new <- births %>% 
    left_join(deaths, by="yb") %>%
    left_join(nmigr %>% rename(yd=ym), by=c("yb","yd","age")) %>% 
    mutate(t = yd + delta) %>% 
    group_by(yb, B, t) %>% 
    summarise(D = sum(D), NM = sum(NM)) %>%
    group_by(yb, B) %>% 
    mutate(Dcum = cumsum(D), NMcum = cumsum(NM)) %>% 
    ungroup() %>% 
    mutate(age = t-yb-delta, pop = B-Dcum+NMcum)
    
  # actual cohorts
  pop_act <- pop  %>% 
    rename(t0 = t, age0 = age, pop0 = pop) %>% 
    left_join(deaths %>% select(-age), by="yb") %>% 
    left_join(nmigr %>% rename(yd=ym) %>% group_by(yd, yb) %>% summarise(NM=sum(NM)),  by=c("yb","yd")) %>% 
    group_by(age0,pop0,t0,yb,yd) %>% 
    summarise(D = sum(D), NM=sum(NM)) %>%
    group_by(age0,pop0,t0,yb) %>% 
    mutate(Dcum = cumsum(D), NMcum = cumsum(NM)) %>% 
    ungroup() %>% 
    mutate(age = age0+yd-t0+delta, pop = pop0-Dcum+NMcum, t=t0+age-age0)
  
  # join
  pop_out <- rbind(
        pop_new %>% select(t, yb, age, pop),
        pop_act %>% select(t, yb, age, pop),
        pop)%>%  
    arrange(t, age)
  
  return(list(pop_out, pop_new, pop_act, deaths, nmigr))
}

# split triangles deaths
split_triangles_D <-function(deaths, births, alpha_input = NULL){
    
    # interval
    delta <- unique(diff(births$yb))
  
    # split triangles
    if(is.null(alpha_input)){
      IMR_obs <- deaths %>% 
        filter(age == 0) %>% 
        inner_join(births, by = c("yd"="yb")) %>% 
        mutate(IMR = D/B)
      alpha_0 <- IMR_obs %>% 
        bind_cols(alpha = sapply(IMR_obs$IMR,function(x) DemoTools::lt_rule_1a0_cd(IMR=x,Sex = "f",region = "w"))) %>% 
        select(-B,-IMR)
      alpha_rest <- deaths %>% filter(age!=0) %>% mutate(alpha=.5)
      alpha <- bind_rows(alpha_0, alpha_rest) %>% arrange(age, yd)
      deaths_splitted <- bind_rows(
        deaths %>% inner_join(alpha) %>% mutate(yb = yd-age-delta, D = D*alpha),
        deaths %>% inner_join(alpha) %>% mutate(yb = yd-age  , D = D*(1-alpha))) %>% 
        arrange(yd, age) %>% 
        select(-alpha) %>% 
        ungroup()
    }
  return(deaths_splitted)
}

# split triangles netmigr
split_triangles_NM <-function(nmigr, alpha_input = .5, delta = 1){
    nmigr_splitted <- bind_rows(
      nmigr %>% mutate(alpha = alpha_input, yb = ym-age-delta, M = M*alpha),
      nmigr %>% mutate(alpha = alpha_input, yb = ym-age  , M = M*(1-alpha))) %>% 
      arrange(ym, age) %>% 
      select(-alpha) %>% 
      ungroup()
    return(nmigr_splitted)
}
  
