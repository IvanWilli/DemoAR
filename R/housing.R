library(devtools)
load_all("C:/Proyectos/DemoTools")

# function to estimate using viviendas

pop_hous_hat <- function(viv_prev, t_prev, viv_n, t_n, porc_ocup, pop_avg){
  viv_n * porc_ocup * pop_avg
}

