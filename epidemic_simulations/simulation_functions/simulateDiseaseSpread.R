#### Function to simulate infection dynamics given values of alpha, beta, epsilon, Tmax, and Coo

simulateDiseaseSpread <- function(alpha = alpha, beta = beta, epsilon = epsilon, Tmax = Tmax, Coo = Coo){
  #require(GillespieSSA)
  source("R_functions/dispersal_kernel_functions.R")
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
          m <- m + normalizedKernel(distance = d, alpha = alpha)
        }
      } # End iInfected loop
      ## Cumulative force of infection for iiSusceptible plant
      lambda[iiSusceptible] <- lambda[iiSusceptible] + beta*m + epsilon
      ## If cumulative lambda (i.e., disease preassure) exceeds Q (i.e., Sellke threshold) plant becomes infected and receives an infection time in range [t - 1, t]
      if(lambda[iiSusceptible] >= Q[iiSusceptible]){
        Inf_times[iiSusceptible] <- runif(1, min = t-1, max = t)
      }
    } # End iiSusceptible loop
  } # End t loop
  ## Order indices of plants based on their infection time
  Inf_indices <- IDs[order(Inf_times)]
  return(list(Inf_times = Inf_times, 
              Inf_indices = Inf_indices))
}


# #### Produce raster with infection times as values
# ## Get raster dimensions from Coo
# ## Coo[,1] = rows or y coord
# ## Coo[,2] = columns or x coord
# CooRows <- Coo[,1]
# CooNrows <- length(unique(CooRows))
# CooYmn <- min(CooRows)
# CooYmx <- max(CooRows)
# CooColumns <- Coo[,2]
# CooNcols <- length(unique(CooColumns))
# CooXmn <- min(CooColumns)
# CooXmx <- max(CooColumns)
# 
# rasterTimes <- raster(nrows = CooNrows, ymn = CooYmn, ymx = CooYmx,
#                       ncols = CooNcols, xmn = CooXmn, xmx = CooXmx,
#                       vals = Inf_times)
# plot(rasterTimes)
