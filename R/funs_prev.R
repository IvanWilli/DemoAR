# conciliation tools


# ajustar población incial - https://timriffe.github.io/DemoTools/reference/basepop_five.html
adjust_base_pop <- function(pop, births, deaths, sex = "M", date = 2010){
  
  # bpe_female <- basepop_single(
  #   #   location = location,
  #   #   refDate = refDate,
  #   #   Females_single = pop_female_single
  #   # )
  #   
  #   # base de al pirámide con nacimientos previos y sobrevivencia
  #   nLx_date <- DemoTools::downloadnLx(nLx=NULL, location = "Argentina",
  #                                      gender = ifelse(sex == "M", "male", "female"), 
  #                                      nLxDatesIn = date)
  #   births_u5 <- births %>% filter(t >= (date-5))
  #   pop %>% filter(t == date) %>% mutate(age)
  #   
  #   # extinct generation method para 100+
  #   deaths_uOAG <- OPAG() # https://timriffe.github.io/DemoTools/reference/OPAG.html
  #   deaths_uOAG <- deaths %>% filter()
}

# get model lt
# get_model_lt<- function(q0_1 = , q0_5){
#   lt <- lt_model_cdun_combin_single(type = "CD_West", Sex = "f",
#                                     q1 = .04651, q5 = .06027,
#                                     indicator_adult = "45q15", value_adult = .241,
#                                     OAnew = 110)
# }

# residual migration. 
  # https://timriffe.github.io/DemoTools/reference/mig_beta.html
  # https://timriffe.github.io/DemoTools/reference/mig_resid.html

# migration patterns by age
# DemoTools::mig_un_families %>% 
#   ggplot(aes(age, prop, color=mig_sign))+
#   geom_line(aes(linetype=factor(family)))+
#   facet_grid(~sex)

# going further pop: 
  # all with same time interval
  # pend: # extinct generations
          # customize ax
          # IMN as time vector
lexis_pop <- function(pop, births, deaths, nmigr = NULL, sex = "b", ax = NULL, IMN = 1/1.98, NewOAG = NULL){
  
  # births <- births_sex
  # deaths <- deaths_sex
  # pop <- pop_sex
  # nmigr <- NULL
  
  births_age_mother <- births
  births <- births %>% group_by(yb) %>% summarise(B = sum(B))
  delta <- unique(diff(births$yb))
  
  # split death triangles: return yb, yd, age, D    
  deaths <- split_triangles_D(deaths, births)
  
  # split smigr triangles: return yb, yd, age, M    
  if(is.null(nmigr)){
    nmigr <- deaths %>% rename(ym=yd, NM=D) %>% mutate(NM=0)
  }else{
    nmigr <- split_triangles_NM(nmigr)  
  }

  # rename pop and filter posterior events
  pop <- pop  %>% mutate(yb = t - age - delta)
  t0 <- unique(pop$t)
  births <- births %>% filter(yb>=t0)
  deaths <- deaths %>% filter(yd>=t0)
  
  # new cohorts 
  pop_new <- births %>% 
    left_join(deaths, by="yb") %>%
    left_join(nmigr %>% rename(yd=ym), by=c("yb","yd","age")) %>% 
    mutate(t = yd + delta) %>% 
    group_by(yb, B, t) %>% 
    summarise(D = sum(D), NM = sum(NM)) %>%
    group_by(yb, B) %>% 
    mutate(Dcum = cumsum(D), NMcum = cumsum(NM)) %>% ungroup() %>% 
    mutate(age = t-yb-delta, pop = B-Dcum+NMcum)
    
  # actual cohorts (take care extinct)
  pop_act <- pop  %>% 
    rename(t0 = t, age0 = age, pop0 = pop) %>% 
    left_join(deaths, by="yb") %>% 
    left_join(nmigr %>% rename(yd = ym) %>% group_by(yd, yb) %>% summarise(NM=sum(NM)),  by=c("yb","yd")) %>% 
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
        pop) %>%
    mutate(age = ifelse(age>NewOAG,NewOAG,age)) %>%
    group_by(t, yb, age) %>% 
    summarise(pop = sum(pop, na.rm=T)) %>% ungroup() %>% 
    mutate(pop = pmax(pop, .5)) %>% # take care of this
    arrange(t, age)
  
  # out
  return(list(pop_out=pop_out, pop_new=pop_new, pop_act=pop_act, 
              deaths=deaths, nmigr=nmigr))
}

# get asfr
get_period_asfr <- function(births_age_mother, pop_f){
  
  # exposures
  t_exposures <-  zoo::rollmean(unique(pop_f$t), k = 2)
  exposures <- interp_tidy(pop_f, age, t, pop, t_exposures) %>% rename(age = x, t=y, exp =value)
  
  # period fertility rates
  asfr <- births_age_mother %>% mutate(t = yb + .5) %>% select(t, age, B) %>%
      right_join(exposures) %>% 
      mutate(asfr = B/exp)
  
  # tfr
  tfr <- asfr %>% group_by(t) %>% summarise(TFR = sum(asfr, na.rm = T))
  
  # pert_asfr
  perc_asfr <- asfr %>% group_by(t) %>% mutate(perc_asfr = asfr/sum(asfr,na.rm=T))
    
  # out
  return(list(asfr=asfr, tfr=tfr, perc_asfr=perc_asfr))
}

# split triangles deaths
split_triangles_D <-function(deaths, births, alpha_input = NULL, max_age = 130){
    
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
      nmigr %>% mutate(alpha = alpha_input, yb = ym-age-delta, NM = NM*alpha),
      nmigr %>% mutate(alpha = alpha_input, yb = ym-age  , NM = NM*(1-alpha))) %>% 
      arrange(ym, age) %>% 
      select(-alpha) %>% 
      ungroup()
    return(nmigr_splitted)
}

# interpolate tidy data
interp_tidy <- function(data, x, y, z, new_y, method = "exponential"){
  x <- enquo(x);y <- enquo(y);z <- enquo(z)
  data_matrix <- data %>% 
      group_by(!!x, !!y) %>% summarise(value = sum(!!z, na.rm=T)) %>% ungroup() %>% 
      pivot_wider(names_from=!!y, values_from=value) %>% as.matrix()
  data_interp <- tibble(x = data_matrix[,1],
                        interp(popmat = data_matrix[,-1],
                             datesIn = as.numeric(colnames(data_matrix)[-1]),
                             datesOut = new_y,
                             method = method, 
                             extrap = T) %>% as.data.frame()) %>%
                      setNames(c("x",as.character(new_y))) %>% 
                      pivot_longer(!x, names_to = "y", values_to = "value") %>%
                      mutate(y = as.numeric(y)) # %>% select(as_label(x) = 1, as_label(y) = 2, as_label(z) = 3)
  return(data_interp)
}

# get asmr
get_period_asmr_lt <- function(deaths, pop, smooth = FALSE, ...){
  
  # exposures
  t_exposures <-  zoo::rollmean(unique(pop$t), k = 2)
  exposures <- interp_tidy(pop, age, t, pop, t_exposures) %>% rename(age = x, t=y, exp =value)
  
  # period mortality rates
  asmr <- deaths %>% 
    group_by(yd, age) %>% summarise(D = sum(D)) %>% ungroup() %>% 
    mutate(t = yd + .5) %>% select(-yd) %>% 
    right_join(exposures) %>% 
    mutate(asmr = D/exp)
  
  # period life tables
  lt <- asmr %>% 
    split(list(.$t)) %>% 
    lapply(., function(X){
      X_t <- unique(X$t)
      X_lt <- DemoTools::lt_single_mx(X$asmr, OAnew = max(X$age), extrapFrom = 90)
      return(X_lt %>% mutate(t = X_t, .before = everything()))
    }) %>% 
    bind_rows()
  
  # out
  out <- list(asmr = asmr, lt = lt)
  return(out)
}

# redistribute unknowns
red_unk <- function(x, v){
  x %>%
    group_by(yb, SEXO) %>%
    mutate(B_star = B + B/sum(B[MEDAD  %in% 10:49], na.rm=T)*sum(B[!MEDAD %in% 10:49], na.rm=T)) %>%
    filter(MEDAD %in% 10:49) %>% ungroup %>%
    select(yb, MEDAD, B = B_star)
}
