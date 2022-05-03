library(devtools)
load_all("C:/Proyectos/DemoTools")

# función para crecimiento vegetativo
# example:
# N0 = 100
# t0 = 2010
# B  = data.frame(B = c(2,3,4,6), t_init = c(2010,2011,2012,2013), t_length = c(1,1,1,.5))
# D  = data.frame(D = c(1,2,2,3), t_init = c(2010,2011,2012,2013), t_length = c(1,1,1,.5))
# cv_total(N0,t0,B,D)
# cv_total(N0,t0,B,D,t_star=2020.4)

cv_total <- function(N0, t0, B, D, t_star=NULL){
  stopifnot(length(B$t_length)==length(D$t_length))
  stopifnot(all(diff(B$t_init)==B$t_length[-nrow(B)]))
  portions <- data.frame(B = pmin(B$t_init+B$t_length - t0, B$t_length),
                         D = pmin(D$t_init+D$t_length - t0, D$t_length))
  N = c(N0 , N0 + cumsum(B$B*portions$B) - cumsum(D$D*portions$D))
  t = c(t0, t0 + cumsum(portions$B)) 
  out <- data.frame(t,N)
  if(!is.null(t_star)){
    out_star <- DemoTools::interp(matrix(out$N,nrow = 2,ncol = length(out$N),byrow = T),
                                         datesIn = out$t, datesOut =  c(out$t,as.numeric(t_star)), 
                                         extrap=T)[1,]
    out <- data.frame(t=c(out$t,as.numeric(t_star)),N=out_star)
    rownames(out) <- NULL
  }
  out
}

# función cv edad, siendoq ue eventos son de la cohorte
cv_edad <- function(N0, t0, B, D, t_star=NULL){
  cv_total(N0, t0, B, D, t_star) 
}