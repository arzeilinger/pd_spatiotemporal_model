#### Function for time-dependent vector density
## Assume damped oscillation form, based phenomenologically on data


vectorDensity <- function(A, lambda_osc, phi, base = 0.0005, time){
  rho <- A*exp(-lambda_osc*time)*sin(phi*pi*time)
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


########################################################################################
#### Exploring seasonal birth rates function of Bilal and Michael 2018
# rho0 <- 1000
# phi <- 0.1
# time <- seq(1,365,by = 1)
# rho <- rho0*(1 + phi*sin(3*pi*time/365))
# 
# #library(ggplot2)
# qplot(time,rho,geom="path", xlab="time", ylab="rho_et")



########################################################################################
#### Function for time-dependent acquisition
#### Should mirror a_m parameter values used in Daugherty and Almeida 2019

## Try a Holling Type IV function

holling4 <- function(a, b, c, time){
  (a*time^2)/(b + c*time + time^2)
}

# a <- 0.1 # long-term y value
# ## x value at peak = -2b/c should coincide with Sept/Oct
# b <- 36*2
# c <- -8*2
# time <- seq(1,12,0.1)

# at <- holling4(a, b, c, time)
# 
# qplot(time,at,geom="path", xlab="time", ylab="at")
