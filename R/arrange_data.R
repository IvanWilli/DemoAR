# pop
library(readxl)
library(tidyverse)
library(foreign)
library(sqldf)
library(xlsx)
library(reshape2)
library(lubridate)
library(readxl)
library(DemoTools)

# pop
pop <- read_xlsx("data/bd_pad_pobl.xlsx", sheet = "bd") %>% 
        mutate(fecha_dec = decimal_date(as.Date(fecha,"%d/%m/%Y"))) %>% 
        filter(round(fecha_dec,1) %in% c(2001.5,2010.5)) %>% 
        select(sexo, edad, valor, fecha_dec) %>% 
        pivot_wider(names_from=fecha_dec, values_from=valor)

pop <- tibble(pop[c("sexo","edad")],
              pop = interp(popmat = as.matrix(pop[,3:4]),datesIn = as.numeric(colnames(pop)[3:4]),
                    datesOut = 2010,method = "exponential")) %>% 
        right_join(pop, by = c("sexo","edad")) %>% 
        rename(`2010`=pop)

# deaths
files <- list.files("data/DEIS/D",full.names = T)
dbD <- map_df(files[20:30], 
       function(i){ # i <- files[25]
         print(i)
          files_ext <- tolower(str_sub(file.path(i),-3))
          if(files_ext=="dbf"){db_i<-read.dbf(i) 
          }else if(files_ext=="csv"){
            sep = ifelse(str_sub(file.path(i),-8,-5) %in% c(2017,2019),";",",")
            db_i<-read.csv(i,sep = sep,stringsAsFactors = F)
            }
          if("PROVRES" %in% names(db_i)) db_i$PROVRE<-db_i$PROVRES
          db_i$anio<-str_sub(file.path(i),-8,-5) %>% as.integer()
          db_i$fo_mm<-as.integer(colsplit(db_i$FECDEF, '/', names =  c('dd','mm', 'aa'))[,2])
          db_i$fo_aa<-as.integer(colsplit(db_i$FECDEF, '/', names =  c('dd','mm', 'aa'))[,3])
          if(is.null(db_i$FECINS)){
            db_i$fr_mm <- db_i$fr_aa <- NA  
          }else{
            db_i$fr_mm<-as.integer(colsplit(db_i$FECINS, '/', names =  c('dd','mm', 'aa'))[,2])
            db_i$fr_aa<-as.integer(colsplit(db_i$FECINS, '/', names =  c('dd','mm', 'aa'))[,3])
          }
          db_i$EDAD[db_i$UNIEDA==2 | db_i$UNIEDA==3 | db_i$UNIEDA==4 | db_i$UNIEDA==8]<-0
          db_i$SEXO[db_i$SEXO==9] <- 3
          db_i$PROVRE<-as.integer(db_i$PROVRE)
          db_i$EDAD<-as.integer(db_i$EDAD)
          db_i$SEXO<-as.integer(db_i$SEXO)
          db_i <- db_i %>% group_by(PROVRE, EDAD, SEXO, fr_mm, fr_aa,  fo_mm, fo_aa, anio) %>% summarise(n = n()) %>% ungroup()
  }
)
write.csv(dbD, "data/D.csv", row.names = F)
dbD %>% group_by(anio) %>% summarise(sum(n))

# births
files <- list.files("data/DEIS/B",full.names = T)[1:29]
dbB <- map_df(files[19:29], 
        function(i){ # i <- files[29]
          print(i)
          files_ext <- tolower(str_sub(file.path(i),-3))
          if(files_ext=="dbf"){db_i<-read.dbf(i,) 
          }else if(files_ext=="csv"){
            sep = ifelse(as.integer(str_sub(file.path(i),-8,-5)) %in% c(2014:2017,2019),";",",")
            db_i<-read.csv(i,sep = sep)
          }
          if("MPRORES" %in% names(db_i)) db_i$PROVRE <-db_i$MPRORES
          if("MPROVRES" %in% names(db_i)) db_i$PROVRE <-db_i$MPROVRES
          if("PROVRES" %in% names(db_i)) db_i$PROVRE <-db_i$PROVRES
          if( "MPRORE" %in% names(db_i)) db_i$PROVRE <-db_i$MPRORE
          if(!"MEDAD" %in% names(db_i)) db_i$MEDAD <-db_i$EDAD
          if( "FECNAC" %in% names(db_i)) db_i$FECOC <-db_i$FECNAC
          if(!"FECOC" %in% names(db_i)) {
            if("MESOC" %in% names(db_i)) db_i$fo_mm<- db_i$MESOC else db_i$fo_mm<- NA
            if("ANO" %in% names(db_i)) db_i$fo_aa<- db_i$ANO else db_i$fo_aa<- NA
          }else{
            db_i$fo_mm<-as.integer(colsplit(db_i$FECOC, '/', names =  c('dd','mm', 'aa'))[,2])
            db_i$fo_aa<-as.integer(colsplit(db_i$FECOC, '/', names =  c('dd','mm', 'aa'))[,3])
          }
          if(is.null(db_i$FECINS)){
            db_i$fr_mm <- db_i$fr_aa <- NA  
          }else{
            db_i$fr_mm<-as.integer(colsplit(db_i$FECINS, '/', names =  c('dd','mm', 'aa'))[,2])
            db_i$fr_aa<-as.integer(colsplit(db_i$FECINS, '/', names =  c('dd','mm', 'aa'))[,3])
          }
          db_i$anio<-str_sub(file.path(i),-8,-5) %>% as.integer()
          db_i$SEXO <- as.integer(db_i$SEXO)
          db_i$SEXO[db_i$SEXO==9] <- 3
          db_i$PROVRE<-as.integer(db_i$PROVRE)
          db_i$MEDAD<-as.integer(db_i$MEDAD)
          db_i <- db_i %>% group_by(PROVRE, MEDAD, SEXO, fr_mm, fr_aa,  fo_mm, fo_aa, anio) %>% summarise(n = n()) %>% ungroup()
        }
)
write.csv(dbB, "data/B.csv", row.names = F)
dbB %>% group_by(anio) %>% summarise(sum(n))

# union 2020
dbD<-read.csv("data/D.csv")
dbD_2020<-read.csv("data/DEIS/IvanDef2020.csv") %>% 
  mutate(fo_mm = as.integer(substr(mesdef,1,2)),
         fo_aa = as.integer(substr(mesdef,4,7)),
         fr_mm = as.integer(substr(mesreg,1,2)),
         fr_aa = as.integer(substr(mesreg,4,7)),
         EDAD = ifelse(uniedad %in% c(2,3,4,8), 0, edad)) %>% 
  rename(anio=`ï..ano`, PROVRE=provres, SEXO=sexo, n=cant) %>% 
  group_by(PROVRE, EDAD, SEXO, fr_mm, fr_aa,  fo_mm, fo_aa, anio) %>% summarise(n = sum(n)) %>% ungroup()
dbD <- bind_rows(dbD %>% filter(anio<2019), dbD_2020) 
dbD %>% group_by(anio) %>% summarise(sum(n))
dbB_2020<-read.csv("data/DEIS/IvanNac2020.csv") %>% 
  mutate(fo_mm = as.integer(substr(mesnac,1,2)),
         fo_aa = as.integer(substr(mesnac,4,7)),
         fr_mm = as.integer(substr(mesreg,1,2)),
         fr_aa = as.integer(substr(mesreg,4,7))) %>% 
  rename(anio=`ï..ano`, PROVRE=provres, SEXO=sexo, MEDAD=medad, n=cant) %>% 
  group_by(PROVRE, MEDAD, SEXO, fr_mm, fr_aa,  fo_mm, fo_aa, anio) %>% summarise(n = sum(n)) %>% ungroup()
dbB<-read.csv("data/B.csv")
dbB <- bind_rows(dbB %>% filter(anio<2019), dbB_2020) 
write.csv(dbB, "data/B.csv", row.names = F)
dbB %>% group_by(anio) %>% summarise(sum(n))

# prepare for export nice
write.xlsx()
dbB %>% group_by(fo_aa) %>% summarise(sum(n))
dbD %>% group_by(fo_aa) %>% summarise(sum(n))
