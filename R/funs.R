# conciliation tools


# redistribute OAG to 100+: observed or with stationary
redistr_OAG <-function(pop, red_from = 80, sex = "Male", e0){
  # method 1: get mlt and red using Lx (preferable using e65 or 45q15)
  
  # method 2: apply observed distribution and smooth
  
  # method 3: splines to omega with lx= .5
  
}

# ecuación de balance total
beq_sex <- function(pop, deaths, births, nmigr = NULL, date_start = unique(pop$t)[1]){
  
  if(is.null(pop$codgeo)){
    pop$codgeo <- NA
    deaths$codgeo <- NA
    births$codgeo <- NA
    nmigr$codgeo <- NA
  }
  if(is.null(nmigr)){
    nmigr <- deaths %>% mutate(NM = 0) %>% select(codgeo, ym = yd, sex, NM)
  }
  
  beq <- births %>% rename(t=yb) %>% 
    left_join(deaths %>% rename(t=yd)) %>% 
    left_join(nmigr %>% rename(t=ym)) %>%
    left_join(pop %>% select(codgeo, sex, pop)) %>% 
    filter(t>=date_start) %>% 
    arrange(codgeo, sex, t) %>% 
    group_by(codgeo, sex) %>% 
    mutate(Dcum = cumsum(D), NMcum = cumsum(NM), Bcum = cumsum(B), 
           new_pop = pop - Dcum + NMcum + Bcum, 
           new_t = t + 1)
  pop_out <- bind_rows(pop, beq %>% select(codgeo, sex, t = new_t, pop = new_pop))
  return(pop_out)
}

# ecuación de balance por edad
beq_by_sex_age <- function(pop, deaths, births, nmigr = NULL, 
                           date_start = unique(pop$t)[1], pop_star = NULL, adjust_init = FALSE){
  
  # pop = pop2010 %>% rename(sex = sexo, age = edad) %>% mutate(age = as.integer(age))
  # deaths = D %>% rename(sex = sexo, age = edad)
  # births = B %>% group_by(sex = sexo, yb) %>% summarise(B = sum(B))
  # nmigr = netmigr_edad_sexo %>% rename(age = edad, NM = nm, ym = t)
  
  results_sex <- list()
  
  for(this_sex in c("V","M")){
    # this_sex = "M"
    
    # pop sex
    pop_sex <- pop %>% filter(sex == this_sex) %>% select(-sex) 
    
    # gratue to single
    if(!is_single(pop_sex$age)){
      pop_sex <- DemoTools::graduate(pop_sex)  
    }
    
    # adjust piramyd - PEND
    if(adjust_init){
      pop_sex <- adjust_base_pop(pop_sex, births, deaths, Px)
    }
    
    # events
    births_sex <- births %>% filter(sex == this_sex, yb>=date_start) %>% select(-sex)
    deaths_sex <- deaths %>% filter(sex == this_sex, yd>=date_start) %>% select(-sex)
    if(!is.null(nmigr)){
      nmigr_sex <- nmigr %>% filter(sex == this_sex, ym>=date_start) %>% select(-sex)
    }else{
        nmigr <- NULL}
    
    # lexis
    pop_out_sex <- lexis_pop(pop_sex, births_sex, deaths_sex, nmigr_sex, sex = this_sex, NewOAG = 100)
    
    # mortality
    asmr_lt_sex <- get_period_asmr_lt(deaths_sex, pop_out_sex$pop_out)
    
    # fertility
    if(this_sex == "M" & !is.null(births$age)) {
      asfr <- get_period_asfr(births, pop_out_sex$pop_out)
    } else{
      asfr <- NULL
      message("no age for asfr")    
    }
    
    # if there is a residual
    if(!is.null(pop_star)){
      show_residual_pop <- show_residual(pop_star, pop_out_sex$pop_out)
      red_residual_pop <- redistr_residual(pop_star, pop_out_sex$pop_out, timing, age_pattern)  
    }
    
    # results
    results_sex[[this_sex]] <- list(pop_out = pop_out_sex, mort_out = asmr_lt_sex, asfr = asfr)
  }
  return(results_sex)
}

# going further pop: # customize ax. IMN as time vector
lexis_pop <- function(pop = NULL, births, deaths, nmigr = NULL, sex = "b",
                      ax = NULL, IMN = 1/1.98, NewOAG = 100, t0 = unique(pop$t)){
  
  # pop = pop2010 %>% rename(sex = sexo, age = edad) %>% mutate(age = as.integer(age)),
  # deaths = D %>% rename(sex = sexo, age = edad),
  # births = B %>% group_by(sex = sexo, yb) %>% summarise(B = sum(B)),
  # nmigr = netmigr_edad_sexo %>% rename(age = edad, NM = nm)
  
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
  deaths <- deaths %>% filter(yd>=t0)
  
  # new cohorts 
  if(!is.null(births)){
    births <- births %>% filter(yb>=t0)
    pop_new <- births %>% 
      left_join(deaths, by="yb") %>%
      left_join(nmigr %>% rename(yd=ym), by=c("yb","yd","age")) %>% 
      mutate(t = yd + delta) %>% 
      group_by(yb, B, t) %>% 
      summarise(D = sum(D), NM = sum(NM)) %>%
      group_by(yb, B) %>% 
      mutate(Dcum = cumsum(D), NMcum = cumsum(NM)) %>% ungroup() %>% 
      mutate(age = t-yb-delta, pop = B-Dcum+NMcum)
  }else{
    pop_new <- pop <- pop_new %>% slice(0)
  }
  
  # actual cohorts (take care extinct)
  if(!is.null(pop)){
    pop <- pop  %>% mutate(yb = t - age - delta)
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
  }else{
    pop_act <- pop <- pop_new %>% slice(0)
  }
  
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

# vase pop 0-10 según mort, fec década pasada - https://timriffe.github.io/DemoTools/reference/basepop_five.html
adjust_base_pop <- function(pop, births, deaths, sex = "M", date = 2010){
  
  refDate <- 1986
  location <- "Brazil"
  pop_female_single <- fertestr::FetchPopWpp2019(location,
                                                 refDate,
                                                 ages = 0:100,
                                                 sex = "female")
  pop_female_counts <- single2abridged(setNames(pop_female_single$pop,
                                                pop_female_single$ages))
  pop_male_single   <- fertestr::FetchPopWpp2019(location,
                                                 refDate,
                                                 ages = 0:100,
                                                 sex = "male")
  pop_male_counts   <- single2abridged(setNames(pop_male_single$pop,
                                                pop_male_single$ages))
  Age <- names2age(pop_male_counts)
  # Automatically downloads the nLx, ASFR, and SRB data
  bpe <- basepop_five(
    location = location,
    refDate = refDate,
    Females_five = pop_female_counts,
    Males_five = pop_male_counts,
    Age = Age
  )
  
  # The counts for the first three age groups have been adjusted:
  bpe$Females_adjusted[1:3]
  pop_female_counts[1:3]
  
  bpe$Males_adjusted[1:3]
  pop_male_counts[1:3]
}

# extinct generation method: D should have yd, age, D
extinct_generation <-function(deaths, births, age_start_extinct = 90, new_OAG = max(deaths$age)){
  # deaths <- D %>% filter(sexo=="V", yd>=2001) %>% select(yd, age=edad, D)
  
  # is enough?
  omega_assumed <- max(deaths$age)
  if(omega_assumed<110)stop("seems not enough years to complete cohort")
  start_year <- min(deaths$yd)
  end_year <- max(deaths$yd)
  development_years <- diff(range(deaths$yd))
  
  # split in case is not
  if(is.null(deaths$yb)){
    deaths <- split_triangles_D(deaths)
  }
  
  # cumulative deaths bu cohort
  new_pop <- deaths %>% 
    mutate(age_start = start_year - yb,
           end_age = age_start + development_years) %>% 
    group_by(yb, age_start) %>% 
    summarise(D = sum(D)) 
    
  # new_pop %>% filter(age_start>=90) %>% summarise(sum(D))
  
  message(paste0(development_years, " years of development. date: end last obs year"))
  return(new_pop %>% 
           ungroup %>% 
           filter(age_start>=age_start_extinct) %>% 
           select(age=age_start, yb, pop=D))
}
 
# split triangles deaths
split_triangles_D <-function(deaths, births = NULL, alpha_input = NULL){
  
  # interval
  delta <- unique(diff(unique(deaths$yd)))
  
  # split triangles
  if(is.null(alpha_input)){
    alpha_rest <- deaths %>% filter(age!=0) %>% mutate(alpha=.5)
    if(min(deaths$age) == 0 & !is.null(births)){
      IMR_obs <- deaths %>% 
        filter(age == 0) %>% 
        inner_join(births, by = c("yd"="yb")) %>% 
        mutate(IMR = D/B)
      alpha_0 <- IMR_obs %>% 
        bind_cols(alpha = sapply(IMR_obs$IMR,function(x) DemoTools::lt_rule_1a0_cd(IMR=x,Sex = "f",region = "w"))) %>% 
        select(-B,-IMR) 
      alpha <- bind_rows(alpha_0, alpha_rest) %>% arrange(age, yd)
    }else{
      message("look not including some age")
      alpha <- alpha_rest
    }
    }
    
    # join deaths
    deaths_splitted <- bind_rows(
      deaths %>% inner_join(alpha) %>% mutate(yb = yd-age-delta, D = D*alpha),
      deaths %>% inner_join(alpha) %>% mutate(yb = yd-age  , D = D*(1-alpha))) %>% 
      arrange(yd, age) %>% 
      select(-alpha) %>% 
      ungroup()
    
  return(deaths_splitted)
}

# survive births
surv_births <- function(births, lt = NULL, deaths = NULL, t_end = NULL){
  
  # births <- B %>% filter(between(yb,2005,2009))%>% rename(age = edad)
  # deaths <- D%>% rename(age = edad)
  # mx <- DemoToolsData::WPP2019_lt %>% filter(LocID == 32, Year == 2008, Sex == "m")
  # lt <- DemoTools::lt_abridged2single(nMx = mx$mx, Age = mx$AgeStart)
  # t_end  <- 2010

  start_year = min(births$yb)
  N_Pb <- NULL
  if(!is.null(lt)){
    lt$nPb = lt$nLx/100000
    N_Pb <- births %>% 
      mutate(Age = t_end - yb - 1) %>% 
      left_join(lt) %>% 
      mutate(N = B * nPb) %>% 
      select(Age, N) %>% 
      filter(!is.na(N)) %>% 
      group_by(Age, yb) %>% 
      summarise(pop = sum(N))%>% 
      select(yb, age = Age, pop)
  }
  
  N_eqb <- NULL
  if(!is.null(deaths)){
    N_eqb <- lexis_pop(pop = NULL, 
                       births = births, 
                       deaths = deaths, 
                       t0 = start_year)$pop_new %>% 
      ungroup %>% 
      filter(t == t_end) %>% 
      select(yb, age, pop)
  }
  
  return(list(N_Pb = N_Pb, N_eqb = N_eqb))
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
                        DemoTools::interp(popmat = data_matrix[,-1] %>% apply(2, as.numeric),
                             datesIn = as.numeric(colnames(data_matrix)[-1]),
                             datesOut = new_y,
                             method = method, 
                             extrap = T,
                             negatives = T) %>% as.data.frame()) %>%
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


to_leslie <- function(Lm, Lf, f, IM = .5, only_sex = "b"){
  
  # Nm <- sample(1:1000, 6)
  # Nf <- sample(1:1000, 6)
  # Lm <- c(99, 65, 33, 22, 11, 5)
  # Lf <- Lm * .8
  # f <- c(.01, .02, .03, .01, .000, .000)
  # IM = .5
  
  stopifnot(all(length(Lm)==length(Lf), length(Lm)==length(f)))
  ages <- length(Lm)
  OAG <- ages-1
  Mf <- Mm <- matrix(0, ages, ages)
  
  Pf <- Lf[-1]/Lf[-ages]
  Pf_OAG <- Lf[ages]/sum(Lf[c(ages-1, ages)])
  Pfb <- Lf[1]/100
  bf <- Pfb * (f + f/2) * (1-IM)
  Mf[1,] <- bf
  diag(Mf[-1,]) <- Pf
  Mf[ages, (ages-1):ages] <- Pf_OAG
  Pm <- Lm[-1]/Lm[-ages]
  Pm_OAG <- Lm[ages]/sum(Lm[c(ages-1, ages)])
  Pmb <- Lm[1]/100
  bm <- Pmb * (f + f/2) * IM
  Bm <- matrix(0, ages, ages)
  Bm[1,] <- bm
  diag(Mm[-1,]) <- Pm
  Mm[ages, (ages-1):ages] <- Pm_OAG
  
  M0 <- matrix(0, ages, ages)
  M <- rbind(cbind(Mf, M0),
             cbind(Bm, Mm))
  
  as.integer(M %*% c(Nf, Nm))
  
  return(M)
}

# do lesliematrix
do_leslie <- function(Sx, Fx){
  ages <- length(Sx)
  leslie <- matrix(0, ages, ages, T)
  leslie[1,] <- Fx
  leslie[row(leslie)-1 == col(leslie)] = Sx[-1]
  leslie
}

# projecting leslie: https://github.com/ubasellini/EDSD2022-population-projections
pop_proj_leslie <- function(N0, Sx, Fx, SMx, IM = 1.04){
  Sx_matrix <- Sx %>% select(t, Sx, Age) %>% pivot_wider(names_from = t,  values_from = Sx) %>% select(-Age)
  Fx_matrix <- Fx %>% pivot_wider(names_from = t, values_from = asfr1) %>% select(-age)
  max_age <- max(Sx$Age)
  N0_vector <- N0 %>% mutate(age = pmin(age, max_age)) %>% group_by(age) %>% summarise(pop = sum(pop)) %>% select(pop)
  P_new <- map2(.x = Sx_matrix, .y = Fx_matrix * IM/(1+IM), ~do_leslie(.x, .y)) %>% 
    accumulate(function(x, y) y %*% x) %>% 
    map_df(~.x %*% N0_vector$pop)
  P_new
}

# projecting tidy
pop_proj_tidy <- function(N0, Sx, Fx, SMx = NULL, IM = 1.04, t_project = 5){
  date_start <- unique(N0$t)
  ages <- sort(unique(N0$Age))
  oag <- max(ages)
  if(is.null(SMx)){
    SMx <- Sx %>% select(t, Age, Sex) %>% distinct() %>% mutate(SMx = 0)
  }
  pop_proj <- list(); births_proj <- list(); deaths_proj <- list(); pop_proj[[1]] <- N0
  for(y in 1:t_project){
      # y = 2
      pop_prev <- pop_proj[[y]]
      # females
      pop_f <- pop_prev %>% filter(Sex == "f") %>%  
              left_join(Sx %>% filter(trunc(t) == (date_start + (y-1))) %>% select(Age, Sex, nLx, Tx)) %>% 
              left_join(Fx %>% filter(trunc(t) == (date_start + (y-1))) %>% mutate(Sex = "f") %>% select(Age, asfr)) %>% 
              left_join(SMx %>%filter(trunc(t) == (date_start + (y-1))) %>% select(Age, Sex, SMx))
      pop_f$Px <- c(pop_f$nLx[1]/100000, 
                  pop_f$nLx[2:oag]/pop_f$nLx[1:(oag-1)],
                  pop_f$Tx[oag+1]/pop_f$Tx[oag])
      pop_f$pop_new <- 0
      pop_f$pop_new[2:oag] <- (pop_f$pop[1:(oag-1)]+pop_f$SMx[1:(oag-1)]/2) * pop_f$Px[2:oag] + pop_f$SMx[2:(oag)]/2
      B_f <- tibble(t0 = date_start + (y-1), t1 = date_start + y,
                    Sex = "f", births = sum((pop_f$pop + pop_f$pop_new + pop_f$SMx/2)/2 * pop_f$asfr) * 1/(1+IM))
      pop_f$pop_new[1] <- B_f$births * pop_f$Px[1] + pop_f$SMx[1]/2
      pop_f$pop_new[oag+1] <- (sum(pop_f$pop[oag:(oag+1)]) + sum(pop_f$SMx[oag:(oag+1)])) * pop_f$Px[oag] + pop_f$SMx[oag+1]/2
      D_f <- tibble(t0 = date_start + (y-1), t1 = date_start + y, Sex = "f", Age = ages, 
                    deaths = c(B_f$births,pop_f$pop[1:(oag-1)], sum(pop_f$pop[oag:(oag+1)])) - pop_f$pop_new)
      # males
      pop_m <- pop_prev %>% filter(Sex == "m") %>%  
              left_join(Sx %>% filter(trunc(t) == (date_start + (y-1))) %>% select(Age, Sex, nLx, Tx)) %>% 
              left_join(SMx %>%filter(trunc(t) == (date_start + (y-1))) %>% select(Age, Sex, SMx))
      pop_m$Px <- c(pop_m$nLx[1]/100000, 
                    pop_m$nLx[2:oag]/pop_m$nLx[1:(oag-1)],
                    pop_m$Tx[oag+1]/pop_m$Tx[oag])
      pop_m$pop_new <- 0
      pop_m$pop_new[2:oag] <- (pop_m$pop[1:(oag-1)]+pop_m$SMx[1:(oag-1)]/2) * pop_m$Px[2:oag] + pop_m$SMx[2:(oag)]/2
      B_m <- tibble(t0 = date_start + (y-1), t1 = date_start + y,
                    Sex = "m", births = sum((pop_f$pop + pop_f$pop_new + pop_f$SMx/2)/2 * pop_f$asfr) * IM/(1+IM))
      pop_m$pop_new[1] <- B_m$births * pop_m$Px[1] + pop_m$SMx[1]/2
      pop_m$pop_new[oag+1] <- (sum(pop_m$pop[oag:(oag+1)]) + sum(pop_m$SMx[oag:(oag+1)])) * pop_m$Px[oag] + pop_m$SMx[oag+1]/2
      D_m <- tibble(t0 = date_start + (y-1), t1 = date_start + y, Sex = "m", Age = ages, 
                    deaths = c(B_m$births,pop_m$pop[1:(oag-1)], sum(pop_m$pop[oag:(oag+1)])) - pop_m$pop_new)
      # get both
      pop <- bind_rows(pop_f %>% mutate(t = t + 1) %>% select(t, Age, Sex, pop = pop_new),
                       pop_m %>% mutate(t = t + 1) %>% select(t, Age, Sex, pop = pop_new))
      deaths <- bind_rows(D_f, D_m) 
      births <- bind_rows(B_f, B_m) 
      pop_proj[[y+1]] <- pop
      deaths_proj[[y]] <- deaths
      births_proj[[y]] <- births
  }
  pop <- bind_rows(pop_proj)
  deaths <- bind_rows(deaths_proj)
  births <- bind_rows(births_proj)
  # out
  return(list(pop=pop, deaths=deaths, births=births))
}


# lee carter function
lc <- function(input, dates_out, jump_off = TRUE, LC_fit_e0 = FALSE, ...){
  dates_in <- sort(unique(as.numeric(input$Date)))
  ages <- sort(unique(input$Age))
  M <- input %>%
    select(Date,nMx,Age) %>% 
    pivot_wider(names_from = Date, values_from = nMx)   %>% 
    select(-Age) %>% as.matrix()
  ndates_in <- ncol(M)
  ax  <- rowSums(log(M))/ndates_in
  M_svd      <- svd(log(M)-ax)
  bx         <- M_svd$u[, 1]/sum(M_svd$u[, 1])
  bx[bx<0] <- 0
  bx         <- bx/sum(bx)
  kt         <- M_svd$d[1] * M_svd$v[, 1] * sum(M_svd$u[, 1]) 
  params <- list(ax = ax, bx = bx, kt = kt)
  # fit k
  if(LC_fit_e0){
    e0_fit <- input$ex[input$Age==0]
    Kt_star <- c()
    for (j in seq_along(dates_in)){
      Kt_star[j] <- optimize(f = interp_lc_lim_kt_min,
                             interval = c(-20, 20),
                             ax = params$ax,
                             bx = params$bx,
                             age = ages,
                             e0_target = e0_fit[j],
                             sex = "b",
                             ...)$minimum
    }
    params$Kt <- Kt_star
  }
  
  kt_diff <- diff(kt)
  summary_kt              <- lm(kt ~ dates_in)
  kt_drift                <- summary_kt$coefficients[2]
  h                       <- dates_out - dates_in[1] # + back_after
  ndates_in <- ncol(M)
  if(all(h<0)){
    if(jump_off){
      ax <- log(M[,1])
    }
    h <- dates_out - min(dates_in)
  }
  if(all(h>0)){
    if(jump_off){
      ax <- log(M[,ndates_in])
    }
    h <- dates_out - max(dates_in)
  }
  kt_forecast <- h * kt_drift
  M_hat <- exp(ax + bx %*% t(kt_forecast)) %>% as.data.frame()
  colnames(M_hat) <- dates_out
  M_hat$Age <- ages
  M_hat <- M_hat %>%
    pivot_longer(cols=-ncol(M_hat), names_to="Date",values_to="nMx")
  M_hat$Sex <- unique(input$Sex)
  out <- M_hat %>% 
    split(list(.$Date, .$Sex), drop = TRUE) %>% 
    lapply(function(X){
      LT <- lt_ambiguous(nMx_or_nqx_or_lx = X$nMx,
                         type = "mx",
                         Age = X$Age,
                         Sex = unique(X$Sex),
                         axmethod = "un",
                         Single = TRUE)
      LT$Sex <- unique(X$Sex)
      LT$Date <- as.numeric(unique(X$Date))
      LT
    }) %>% 
    do.call(rbind,.) 
  return(out)
}

# resultado de api wpp
api_wpp_results <- function(indicator = 2, year_start = 2005, year_end = 2010, location = 32){
  # Declares the base url for calling the API
  base_url <- "https://population.un.org/dataportalapi/api/v1"
  # tomado de https://population.un.org/dataportal/about/dataapi#r-example3-returning-a-single-indicator-for-a-single-geographical-area
  target <- paste0(base_url,"/data/indicators/",indicator,"/locations/",location, "/start/",year_start,"/end/",year_end)
  # Call the API
  response <- jsonlite::fromJSON(target)
  # Get the first page of data
  df <- response$data
  # Get the other pages with data
  while (!is.null(response$nextPage)){
    response <- jsonlite::fromJSON(response$nextPage)
    df <- rbind(df, response$data)
  }
  return(df)
}

# indicadores disponibles wpp
api_wpp_indicadores_disponibles <- function(x){
  # Declares the base url for calling the API
  base_url <- "https://population.un.org/dataportalapi/api/v1"
  # Creates the target URL, indicators, in this instance
  target <- paste0(base_url, "/indicators/")
  # Get the response, which includes data as well as information on pagination and number of records
  response <- jsonlite::fromJSON(target)
  # Get the first page of data
  df <- response$data
  # Loop until there are new pages with data
  while (!is.null(response$nextPage)){
    #call the API for the next page
    response <- jsonlite::fromJSON(response$nextPage)
    #add the data of the new page to the data.frame with the data of the precious pages
    df <- rbind(df, response$data)
  }
  return(df)
}




