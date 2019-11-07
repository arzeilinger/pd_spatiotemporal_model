#### Seasonal PD model
## Based on model of Parry et al. (2014) PNAS

my.packages <- c("ggplot2", "raster")
lapply(my.packages, require, character.only = TRUE)

#### Source additional functions
## Vector density function
source("R_functions/seasonal_vector_density_functions.R")
## Dispersal kernel functions, also includes makeCoordinates() function
source("R_functions/dispersal_kernel_functions.R")


#### Parameter values
## In this version, epsilon and beta are functions of other parameters, define non-varying parameters here
alpha <- 10 # Dispersal parameter
eta <- 0.5 # Inoculation rate (LAMBDA in Parry et al.)
kappa_e <- 0.6 # Proportion of external vectors infectious
a <- 0.6 # Acquisition rate
muv <- 0.2 # In-field loss rate of vectors (due to both death and emigration)

## Number of time steps
## Assume each time step is ~2 weeks
Tmax <- 100

#### Plant spatial coordinates
## Total number of plants (numPlants) = nrow(Coo) = nrc*nrc
nrc <- 50
Coo <- makeCoordinates(nrc)

## Calculating distance to nearest field border for each plant
borderCoo <- rep(c(1,nrc), each = 2)
borderD <- apply(Coo, 1, function(x) min(sqrt((x - borderCoo)^2))+1)
## Add 1 to the distance; plants on the border are then 1 space from the border and all other plants are shifted accordingly
## Dispersal distances from border
## When d = 0, normalizedKernel = Inf
epsilonK <- chooseKernel(d = borderD, alpha = alpha, kernelFunc = "normalized exponential")
qplot(borderD,m_epsilon,geom="path", xlab="distance", ylab="K(D,alpha)")


#### Setting up Inf_times and Inf_indices
## Assumes already created a grid of hosts
## Calculate total number of plants and vector of plant IDs
numPlants <- nrow(Coo)
IDs <- 1:numPlants
## Start with infection times all at Tmax
Inf_times <- rep(Tmax, numPlants)
## Vector of lambda values for each plant
lambda <- rep(0, numPlants)
## Vector of Sellke thresholds
Q <- rexp(numPlants, rate = 1)
## Nested for loops of t, Susceptible plants, and Infected plants
for(t in 1:Tmax){
  ## Each time point, update vectors of infecteds, susceptibles, and their respective indices
  Infected_ind <- which(Inf_times < Tmax)
  #Infected_times <- Inf_times[Infected_ind]
  Susceptible_ind <- which(Inf_times == Tmax)
  #Susceptible_times <- Inf_times[Susceptible_ind]
  ## Summaries of infections and susceptibles
  numInfections <- length(Infected_ind)
  numSusceptibles <- length(Susceptible_ind)
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
        m <- m + chooseKernel(distance = d, alpha = alpha, kernelFunc = "normalized exponential")
      }
    } # End iInfected loop
    ## Calculate time-dependent epsilon
    rho_et <- vectorDensity(A = -1.2, lambda_osc = 0.2, phi = 0.4, base = 0.05, time = t) # Seasonal vector density as damped sine wave
    epsilonti <- eta*kappa_e*rho_et*epsilonK[iiSusceptible]
    ## Calculate time-dependent beta
    Infecteds_tm1 <- sum(Inf_times < t) # Number of infected plants at t-1
    kappa_t <- a*Infecteds_tm1/numPlants # Time-dependent proportion vectors infectious
    rho_tm1 <- (1-muv)*vectorDensity(A = -1.2, lambda_osc = 0.2, phi = 0.4, base = 0.05, time = t - 1)
    beta_t <- eta*kappa_t*rho_tm1
    ## Cumulative force of infection for iiSusceptible plant
    lambda[iiSusceptible] <- lambda[iiSusceptible] + beta_t*m + epsilonti
    ## If cumulative lambda (i.e., disease preassure) exceeds Q (i.e., Sellke threshold) plant becomes infected and receives an infection time in range [t - 1, t]
    if(lambda[iiSusceptible] >= Q[iiSusceptible]){
      Inf_times[iiSusceptible] <- runif(1, min = t-1, max = t)
    }
  } # End iiSusceptible loop
} # End t loop
## Order indices of plants based on their infection time
Inf_indices <- IDs[order(Inf_times)]


#### Produce raster with infection times as values
## Get raster dimensions from Coo
## Coo[,1] = rows or y coord
## Coo[,2] = columns or x coord
CooRows <- Coo[,1]
CooNrows <- length(unique(CooRows))
CooYmn <- min(CooRows)
CooYmx <- max(CooRows)
CooColumns <- Coo[,2]
CooNcols <- length(unique(CooColumns))
CooXmn <- min(CooColumns)
CooXmx <- max(CooColumns)

rasterTimes <- raster(nrows = CooNrows, ymn = CooYmn, ymx = CooYmx,
                      ncols = CooNcols, xmn = CooXmn, xmx = CooXmx,
                      vals = Inf_times)
plot(rasterTimes)

