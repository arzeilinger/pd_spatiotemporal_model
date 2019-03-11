#### FUNCTION TO SIMULATE EPIDEMIC PROCESS
## Result is a simulated data set from infection dynamices
## According to Adrakey's force of infection equation

rm(list = ls())

#### Kernel functions
## Simple exponential kernel
kernel <- function(distance, alpha) {
  exp(-(distance / alpha))
}

## Normalized exponential kernel
## Derived in Neri et al. 2014
## Used in Adrakey et al. 2017
normalizedKernel <- function(distance, alpha) {
  (1/(2*pi*distance*alpha))*exp(-(distance / alpha))
}


#### Initial values setup
Tmax <- 100
## Parameter values used in Adrakey et al. 2017 simulations
## Values are in units day-1 km2
alpha <- 0.08 # km^2
beta <- 7e-6 # day^-1 km^2
epsilon <- 5e-5 # day^-1
grid <- 1:10
Coo <- matrix(c(rep(grid, 10), rep(grid, each = 10)), ncol = 2)
numPlants <- nrow(Coo)
IDs <- 1:numPlants
head(Coo)

## Setting up Inf_times and Inf_indices
## Start with infection times all at Tmax
Inf_times <- rep(Tmax, numPlants)
## For now, start with 1 random plant infected
initial_infected_ind <- sample(1:nrow(Coo), 1)
## Update Inf_times with initial infected
Inf_times[initial_infected_ind] <- runif(1)
## For now, starting with just 1 infected plant from t = 0 to t = 1. Should make "initially infected" a function of epsilon
Inf_indices <- IDs[order(Inf_times)]


for(t in 1:Tmax){
  ## Each time point, update vectors of infecteds, susceptibles, and their respective indices
  Inf_indices <- IDs[order(Inf_times)]
  Infected_ind <- which(Inf_times < Tmax)
  Infected_times <- Inf_times[Infected_ind]
  Susceptible_ind <- which(Inf_times == Tmax)
  Susceptible_times <- Inf_times[Susceptible_ind]
  ## Summaries of infections and susceptibles
  numInfections <- length(Infected_times)
  numSusceptibles <- length(Susceptible_times)
  ## Loop over susceptible plants
  for(iiSusceptible in Susceptible_ind){
    infTimeSusceptible <- Inf_times[iiSusceptible]
    CooSusceptible <- Coo[iiSusceptible,1:2]
    ## (Re)set dispersal kernel summation to 0
    m <- 0
    ## Loop over infected plants to get dispersal kernel from each Infected to a given Susceptible
    for(iInfected in Infected_ind) {
      infTimeInfected <- Inf_times[iInfected]
      if(infTimeSusceptible > infTimeInfected) { ## This avoids Target plant being same as Source
        ## Calculate d for each Infected plant
        d <- sum(sqrt((Coo[iInfected,1:2] - CooSusceptible)^2))
        ## Calculate the dispersal kernel based on each d and then sum over all kernels
        ## Using Neri et al. 2014 normalized exponential kernel
        ## When using simple exponential kernel, m becomes too large, lambdaii > 1, and rbinom fails with error
        m <- m + normalizedKernel(distance = d, alpha = alpha)
      }
    }
    ## Force of infection for iiSusceptible plant
    lambdaii <- beta*m + epsilon
    ## Infection status of iiSusceptible based on lambda
    infectionStatusii <- rbinom(1, 1, lambdaii) 
    ## If iiSusceptible plant is infected, generate a random infection time between t and t+1 time points
    if(infectionStatusii == 1) {
      Inf_times[iiSusceptible] <- runif(1, t, t+1)
    }
  }
}

