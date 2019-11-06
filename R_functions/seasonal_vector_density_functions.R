#### Function for time-dependent vector density
## Assume damped oscillation form, based phenomenologically on data


vectorDensity <- function(A, lambda_osc, phi, base = 0.0005, time){
  rho <- A*exp(-lambda_osc*t)*sin(phi*pi*t)
  ## set negative numbers to zero
  rho <- ifelse(rho < base, base, rho)
  return(rho)
}

#### Parameter values that qualitatively replicate NC PD dynamics:
# A = -1.2
# lambda_osc = 0.2
# phi = 0.4
# t <- seq(0,24,1)
# 
# testrho <- vectorDensity(A = -1.2, lambda_osc = 0.05, phi = 0.2, base = 0.05, t)
# 
# #library(ggplot2)
# qplot(t,testrho,geom="path", xlab="time", ylab="rho_et")

