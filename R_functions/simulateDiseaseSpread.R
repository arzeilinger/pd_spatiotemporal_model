#### Function to simulate infection dynamics given values of alpha, beta, epsilon, Tmax, and Coo

simulateDiseaseSpread <- function(alpha = alpha, beta = beta, epsilon = epsilon, Tmax = Tmax, Coo = Coo){
  require(GillespieSSA)
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
    #Susceptible_times <- Inf_times[Susceptible_ind]
    ## Summaries of infections and susceptibles
    numInfections <- length(Infected_times)
    numSusceptibles <- length(Susceptible_ind)
    ## lambda calculated only for susceptible plants
    lambda <- rep(0, numSusceptibles)
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
      lambda[which(Susceptible_ind == iiSusceptible)] <- beta*m + epsilon
      # ## From Adrakey paper:
      # ## P(ii infected in [t, t+dt]) = lambdaii*dt + o(dt)
      # ## Where dt = period between discrete time points t and t+1
      # ## Assuming that infection status of iiSusceptible is a Bernoulli trial with prob = lambdaii
      # ## lambdaii only stays bounded by [0,1] if normalizedKernel is used
      # infectionStatusii <- rbinom(1, 1, lambdaii) 
      # ## If iiSusceptible plant is infected, generate a random infection time between t and t+1 time points
      # if(infectionStatusii == 1) {
      #   Inf_times[iiSusceptible] <- runif(1, t, t+1)
      # }
    }
    ## Implementing Gillespie SSA within simulation
    params <- lambda
    names(params) <- paste("lambda", 1:numSusceptibles, sep = "")
    a <- names(params)
    x0 <- rep(c(1,0), numSusceptibles)
    names(x0) <- c(paste(c("S", "I"), floor(seq(1, (numSusceptibles+0.5), 0.5)), sep = ""))
    nu <- matrix(0, nrow = length(x0), ncol = length(lambda))
    for(i in 1:ncol(nu)){
      Sind <- paste("S", i, sep = "")
      Iind <- paste("I", i, sep = "")
      nu[which(names(x0) == Sind), i] <- -1
      nu[which(names(x0) == Iind), i] <- 1
    }
    ssaOut <- ssa(x0 = x0,
                  a = a,
                  nu = nu,
                  parms = params,
                  tf = 1,
                  method = ssa.d())
    ## Need to figure out how to specify the time steps
    ssaData <- ssaOut$data
    ssaData[,1:30]
    ssaTimes <- as.numeric(ssaData[,1])
    goodTimes <- ssaTimes[which(ssaTimes <= 1)] + (t - 1)
    ssaData <- ssaData[which(ssaTimes < 1),]
    if(any(ssaTimes > 0 & ssaTimes <=1)){
      Icols <- ssaData[,grep("I", attr(ssaData, "dimnames")[[2]])]
      for(i in 1:ncol(Icols)){
        iName <- attr(Icols, "dimnames")[[2]][i]
        iInd <- as.numeric(gsub("[^0-9]", "", iName))
        isum <- sum(Icols[,i])
        if(isum > 0){
          Inf_times[iInd] <- goodTimes[min(which(Icols[,i] == 1))]
        }
      }
    }
    (infectedNames <- attr(Icols, "dimnames")[[2]][which(colSums(Icols) >= 1)]) 
    # infectedInd <- as.numeric(gsub("[^0-9]", "", infectedNames))
    # Inf_times_ssa <- 
  }
  return(list(Inf_times = Inf_times, 
              Inf_indices = Inf_indices))
}
