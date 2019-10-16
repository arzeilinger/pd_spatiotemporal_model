#### Function to simulate infection dynamics given values of alpha, beta, epsilon, Tmax, and Coo

simulateDiseaseSpread <- function(alpha = alpha, beta = beta, epsilon = epsilon, Tmax = Tmax, Coo = Coo){
  require(GillespieSSA)
  source("R_functions/dispersal_kernel_functions.R")
  #### Setting up Inf_times and Inf_indices
  ## Assumes already created a grid of hosts
  ## Calculate total number of plants and vector of plant IDs
  numPlants <- nrow(Coo)
  IDs <- 1:numPlants
  ## Start with infection times all at Tmax
  Inf_times <- rep(Tmax, numPlants)
  ## Initialize Inf_indices
  Inf_indices <- IDs[order(Inf_times)]
  ## Nested for loops of t, Susceptible plants, and Infected plants
  for(t in 1:Tmax){
    ## Each time point, update vectors of infecteds, susceptibles, and their respective indices
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
      } # End iInfected loop
      ## Force of infection for iiSusceptible plant
      lambda[which(Susceptible_ind == iiSusceptible)] <- beta*m + epsilon
    } # End iiSusceptible loop
    ############################################################################################################
    #### Implementing Gillespie SSA within simulation
    ## Define params, a, and x0 objects for ssa()
    ## Using Inf_indices (precisely, Susceptible_ind) to keep track of the plants
    params <- lambda
    names(params) <- paste("lambda", Susceptible_ind, sep = "")
    a <- names(params)
    x0 <- rep(c(1,0), numSusceptibles)
    names(x0) <- c(paste(c("S", "I"), rep(Susceptible_ind, each = 2), sep = ""))
    ## Define transition matrix, nu
    ## nrow = number of states (S and I for each plant)
    ## ncol = number of rates (lambda for each plant)
    nu <- matrix(0, nrow = length(x0), ncol = length(lambda))
    for(i in 1:ncol(nu)){
      Sind <- paste("S", i, sep = "")
      Iind <- paste("I", i, sep = "")
      nu[which(names(x0) == Sind), i] <- -1
      nu[which(names(x0) == Iind), i] <- 1
    }
    ## Run ssa() function with direct method
    ssaOut <- ssa(x0 = x0,
                  a = a,
                  nu = nu,
                  parms = params,
                  tf = 1,
                  method = ssa.d())
    #### Extract infection states and infection times from output of ssa()
    ## TO DO: Figure out why time steps exceed tf (final time) and how to eliminate these excess time steps
    ssaData <- ssaOut$data
    #print(ssaData[,1:9]) # Check data set, for debugging
    ssaTimes <- as.numeric(ssaData[,1]) # Extract time steps
    goodTimes <- ssaTimes[which(ssaTimes <= 1)] # Filter out time steps that exceed tf = 1
    ssaData <- ssaData[which(ssaTimes <= 1),]
    ## if() statement to catch situations where no infections (i.e., transitions) at t < tf occurred
    if(any(ssaTimes > 0 & ssaTimes <= 1)){
      ## Only the infected state columns
      Icols <- ssaData[,grep("I", attr(ssaData, "dimnames")[[2]])]
      ## For loop to update Inf_times
      for(i in 1:ncol(Icols)){
        iName <- attr(Icols, "dimnames")[[2]][i]
        iInd <- as.numeric(gsub("[^0-9]", "", iName))
        isum <- sum(Icols[,i])
        if(isum > 0){
          Inf_times[iInd] <- goodTimes[min(which(Icols[,i] == 1))] + (t - 1)
        } # End 2nd if() statement
      } # End i loop
    } # End 1st if() statement within Gillespie SSA code
    ## Sort infection indices according to updated infection times
    Inf_indices <- IDs[order(Inf_times)]
    ## To look at which plants became infected, for debugging purposes
    #infectedNames <- attr(Icols, "dimnames")[[2]][which(colSums(Icols) >= 1)]
    #print(infectedNames)
  } # End t loop
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
