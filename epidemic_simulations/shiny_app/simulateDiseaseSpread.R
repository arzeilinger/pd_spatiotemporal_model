#### Function to simulate infection dynamics given values of alpha, beta, epsilon, Tmax, and Coo

simulateDiseaseSpread <- function(alpha = alpha, beta = beta, epsilon = epsilon, Tmax = Tmax, Coo = Coo){
  #### Setting up Inf_times and Inf_indices
  ## Assumes already created a grid of hosts
  ## Calculate total number of plants and vector of plant IDs
  numPlants <- nrow(Coo)
  IDs <- 1:numPlants
  ## Start with infection times all at Tmax
  Inf_times <- rep(Tmax, numPlants)
  ## For now, start with 1 random plant infected prior to t = 1
  ## Should make initial infection dependent on epsilon?
  #initial_infected_ind <- sample(1:nrow(Coo), 1)
  ## Update Inf_times and Inf_indiceswith initial infected
  #Inf_times[initial_infected_ind] <- runif(1)
  Inf_indices <- IDs[order(Inf_times)]
  ## Nested for loops of t, Susceptible plants, and Infected plants
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
      ## From Adrakey paper:
      ## P(ii infected in [t, t+dt]) = lambdaii*dt + o(dt)
      ## Where dt = period between discrete time points t and t+1
      ## Assuming that infection status of iiSusceptible is a Bernoulli trial with prob = lambdaii
      ## lambdaii only stays bounded by [0,1] if normalizedKernel is used
      infectionStatusii <- rbinom(1, 1, lambdaii) 
      ## If iiSusceptible plant is infected, generate a random infection time between t and t+1 time points
      if(infectionStatusii == 1) {
        Inf_times[iiSusceptible] <- runif(1, t, t+1)
      }
    }
  }
  return(list(Inf_times = Inf_times, 
              Inf_indices = Inf_indices))
}
