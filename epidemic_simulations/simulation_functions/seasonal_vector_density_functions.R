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



##########################################################################################
#### Winter recovery function
## Functional form based on Feil et al. (2003) and same as used by Gruber and Daugherty (2013)
hostRecovery <- function(time, b = -0.045){
  ## The curve can be shifted by modifying the value of parameter b
  ## More negative values of b shift curve to the left, less negative values shift to the right.
  ## Default value is estimate from Feil et al. 2003
  (1 - exp(b*time))^135.9
}

# timestep <- 1:365
# timestep/7
# test <- hostRecovery(timestep)
# qplot(x = timestep, y = test, geom = "path")


########################################################################################
#### Function for time-dependent acquisition
#### Should mirror a_m parameter values used in Daugherty and Almeida 2019

## Try a Holling Type IV function

# holling4 <- function(a, b, c, time){
#   (a*time^2)/(b + c*time + time^2)
# }

# a <- 0.1 # long-term y value
# ## x value at peak = -2b/c should coincide with Sept/Oct
# b <- 36*2
# c <- -8*2
# time <- seq(1,12,0.1)

# at <- holling4(a, b, c, time)
# 
# qplot(time,at,geom="path", xlab="time", ylab="at")


########################################################################################
#### Function to process final Exposure times (Exp_times) into cumulative (latent) infections over time
processInfectionTimes <- function(Exp_times, numYears, Tmax, weekVector){
  #### Inputs
  ## Exp_times = a matrix (numPlants x numYears) of the time that each plant entered the Exposed compartment in each year
  ## numYears = scalar number of years
  ## Tmax = scalar number of time steps (weeks) in each year
  ## weekVector = a vector with each consecutive time step (weeks)
  #### Returns a data.frame with columns of time step, the number of new infections at that time, and the cumulative number of infections
  finalInfectionTimes <- ceiling(Exp_times)
  infectionTimesList <- vector("list", numYears)
  for(y in 1:numYears){
    times.y <- finalInfectionTimes[,y]
    times.y <- times.y[times.y < Tmax]
    infCounts <- rle(sort(times.y))
    infectionTimesList[[y]] <- data.frame(timeStep = infCounts$values+(Tmax*(y-1)), 
                                          numInfections = infCounts$lengths,
                                          cumulInfections = cumsum(infCounts$lengths))
    ## Bug: 
    ## If there are no chronic infections, the number of cumulative infections doesn't go to zero at the start of the next year
    ## instead, cumulative infections go to zero at the time of the next infection.
    ## if() statement to fix the bug:
    if(all(Exp_times[,y] != 0)){
      infectionTimesList[[y]] <- rbind(c((Tmax*(y-1)), 0, 0), infectionTimesList[[y]])
    }
  }
  ## Compile infection (exposure) times into a data frame
  infTimesProcessed <- infectionTimesList %>% rbindlist() %>% as.data.frame() %>%
    ## Expand data set so that each time step is included explicitly
    right_join(., weekVector, by = "timeStep") %>%
    fill(cumulInfections) %>% # Replace NAs with previous non-NA value
    ## Replace NAs with 0
    ## Can replace NAs in cumulInfections with 0 now because these are just at the start of the time series
    replace_na(list(numInfections = 0, cumulInfections = 0))
  return(infTimesProcessed)
}


##################################################################################################################
#### Vector overwintering infectivity decay function

VOdecay <- function(kappa_b, a = 4, t){
  kappa_b*exp(-t/a)
}

# kappa_b <- 0.6
# t <- 1:52
# test <- VOdecay(kappa_b, t = t)
# qplot(t/4,test,geom="path", xlab="time", ylab="kappa_beta")



##################################################################################################################
#### Function relating seasonal time of infection and latent/incubation periods

transitionPeriod <- function(a = 8, b = 26, c = 0.02, t){
  ## Use a reparameterized quadratic function: a + bt + ct^2
  ## At minimum, t = b, y = a
  a + c*(t-b)^2
}

# t <- 1:52
# test <- transitionPeriod(a = 4, b = 26, c = 0.02, t = t)
# qplot(t/4,test,geom="path", xlab="Month", ylab="Transition period")
