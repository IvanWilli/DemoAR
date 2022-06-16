# arg case

library(tidyverse)
library(DemoTools)
library(readxl)
library(lubridate)
source("R/funs.R")
lapply(list.files(path = "R/MortpakAbacusLegacyMLT/R", pattern = ".R"), source)

  # load data
deaths <- read.csv("data/D.csv") %>% filter(PROVRE==10) %>%
  mutate(sexo = case_when(SEXO == 1 ~ "V", SEXO == 2 ~ "M", TRUE ~ "U")) %>% 
  group_by(anio, sexo, EDAD) %>% summarise(D = sum(n)) %>% rename(yd = anio, age=EDAD) %>% ungroup()
births <- read.csv("data/B.csv") %>% filter(PROVRE==10) %>%
  mutate(sexo = case_when(SEXO == 1 ~ "V", SEXO == 2 ~ "M", T ~ "U")) %>% 
  group_by(anio, sexo, MEDAD) %>% summarise(B = sum(n)) %>% rename(yb=anio, age=MEDAD) %>% ungroup()
pop <- read_xlsx("data/bd_pad_pobl.xlsx", sheet = "bd") %>%
  mutate(fecha_dec = decimal_date(as.Date(fecha,"%d/%m/%Y"))) %>%
  filter(fuente == "AD35", codgeo==10, round(fecha_dec,1) %in% c(2001.5,2010.5)) %>%
  select(sexo, edad, geo, pop = valor, t = fecha_dec)

# redistr unkn
  # births <- births %>% red_unk(sexo) %>% red_unk(medad)
  # deaths <- deaths %>% red_unk(sexo) %>% red_unk(edad)

# get results
results_sex <- list()
for(sex in c("V","M")){
  # sex = "M"
  # step in 1st jan 2010
  pop_sex <- pop %>% filter(sexo == sex) %>% interp_tidy(., edad, t, pop, 2010) %>% rename(age = x, t = y, pop = value)
  # gratue to single
  pop_sex <- DemoTools::graduate()
  
  # adjust piramyd - PEND
  pop_sex <- adjust_base_pop(pop_sex, births, deaths, Px)
  # events
  births_sex <- births %>% filter(sexo == sex) %>% select(-sexo)
  deaths_sex <- deaths %>% filter(sexo == sex) %>% select(-sexo)
  # lexis
  pop_out_sex <- lexis_pop(pop_sex, births_sex, deaths_sex)
  # mortality
  asmr_lt_sex <- get_period_asmr_lt(deaths_sex, pop_out_sex$pop_out)
  # fertility
  if(sex == "M") {asfr <- get_period_asfr(births, pop_out_sex$pop_out)} else{asfr <- NULL}
  # residual
  show_residual_pop <- show_residual(pop_star, pop_out_sex$pop_out)
  red_residual_pop <- redistr_residual(pop_star, pop_out_sex$pop_out, timing, age_pattern)
  # results
  results_sex[[sex]] <- list(pop_out = pop_out_sex, mort_out = asmr_lt_sex, asfr = asfr)
}
names(results)

# both sex results
