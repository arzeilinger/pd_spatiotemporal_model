#### Function to run seasonal PD model
## Based on model of Parry et al. (2014) PNAS

seasonalPDmodel <- function(parameterList, nrc = nrc, Tmax = Tmax, numYears = numYears, numPlantsRogued = numPlantsRogued, vectorOverwintering = TRUE){
  #### parameterList is a named list specifying the values of the following parameters:
  ## alpha = Dispersal parameter
  ## eta = Inoculation rate (LAMBDA in Parry et al.)
  ## kappa_e0 = Baseline proportion of external vectors infectious
  ## aI = Acquisition rate from Infectious hosts
  ## aD = Acquisition rate from Diseased hosts; aD should be >= aI because of Xylella population growth
  ## muv In-field loss rate of vectors (due to both death and emigration)
  ### Parameter values for vectorDensity() damped sine wave
  ## A = Initial amplitude
  ## lambda_osc = Dampening
  ## phi = Angular frequency
  ## base = Baseline vector density
  #### Vector overwintering
  ## muv_e = Loss rate of vectors overwintering
  #### Winter recovery
  ## b = shifting recovery probability curve
  ##########################################################################################################
  #### Roguing
  ## numPlantsRogued = total number of plants to rogue every year
  ## nrc is the number of plants along a row and column; creates a square field of plants with dimentions nrc x nrc
  ## Tmax = number of time steps within a season/year; currently parameterized for time step at the 1 week scale
  ## numYears = number of years
  ## vectorOverwintering is a logical flag; when TRUE, assumes that initial vector densities and infectivity are....
  ## ... are dependent on densities and infectivity at the end of the previous year.
  ## Load parameters from list to global environment
  list2env(parameterList, globalenv())
  #########################################################################################################
  Coo <- makeCoordinates(nrc)
  ## Calculating distance to nearest field border for each plant
  # borderCoo <- rep(c(1,nrc), each = 2)
  # borderD <- apply(Coo, 1, function(x) min(sqrt((x - borderCoo)^2))+1)
  ## Distance from the left border is just be the 1st coordinate
  borderD <- Coo[,1]
  ## Dispersal distances from border
  epsilonK <- chooseKernel(d = borderD, alpha = alpha, kernelFunc = "normalized exponential")
  #### Setting up empty vectors and matrices to track things
  ## Calculate total number of plants and vector of plant IDs
  numPlants <- nrow(Coo)
  IDs <- 1:numPlants
  #### Set up vectors to track state of each plants
  ## Exposure times should be a matrix: rows = plants, columns = years
  Exp_times <- matrix(Tmax, nrow = numPlants, ncol = numYears)
  ## Vectors of rho_epsilon(t), rho_beta(t), and beta(t) values
  rho_etVec <- rho_btVec <- rep(0, Tmax)
  ## Matrix of epsilon_i(t) values
  epsilonMatrix <- betaMatrix <- matrix(0, nrow = numPlants, ncol = Tmax)
  ## Matrix of total kappa(t) values (kappa_b + kappa_e)
  kappaMatrix <- matrix(0, nrow = Tmax, ncol = numYears)
  #### Start year loop
  for(y in 1:numYears){
    print(paste("On year ", y, sep = ""))
    ####################################################################################################
    #### Stochastic transitions among host compartments
    ## (Re)set Infection times, Disease times (i.e., time of symptom onset)
    Inf_times <- Disease_times <- rep(Tmax, numPlants)
    ## (Re)set lambda values for each plant
    lambda <- rep(0, numPlants)
    ## (Re)set Sellke thresholds for infection
    Q <- rexp(numPlants, rate = 1)
    #### (re)set stochastic latent and cryptic periods
    ## Parry et al. (2014) assumed cryptic period (I -> D) is gamma distributed
    ## Here, we make a similar assumption for latent period (E -> I)
    ## Assume that each period has an average of 4 weeks
    ## a*s (mean) = 4; a:s (ratio) = 4:1; this ratio produces a humped distribution
    ## Parameterization assumes that time step = 1 week, will need to be adjusted if time step changes
    ## Durations of stays in Exposed compartment
    Et <- rgamma(numPlants, shape = 4, scale = 1)
    ## Durations of stays in Infectious compartment
    It <- rgamma(numPlants, shape = 4, scale = 1)
    ## In-field vector infectivity from previous year
    kappa_b_ym1 <- kappaMatrix[Tmax, y-1]
    #####################################################################################################
    #### Nested for loops of t, Susceptible plants, and Infected plants
    for(t in 1:Tmax){
      ## Each time point, update vectors of infecteds, susceptibles, and their respective indices
      Infected_ind <- which(Inf_times < Tmax)
      Susceptible_ind <- which(Exp_times[,y] == Tmax)
      ## Summaries of infections and susceptibles
      numInfections <- length(Infected_ind)
      numSusceptibles <- length(Susceptible_ind)
      #######################################################################################################
      ## Loop over susceptible plants
      for(iiSusceptible in Susceptible_ind){
        infTimeSusceptible <- Exp_times[iiSusceptible, y]
        CooSusceptible <- Coo[iiSusceptible,1:2]
        ## (Re)set dispersal kernel summation to 0
        m <- 0
        #### Loop over infected plants to get dispersal kernel from each Infected to a given Susceptible
        for(iInfected in Infected_ind) {
          infTimeInfected <- Inf_times[iInfected]
          if(infTimeSusceptible > infTimeInfected) { ## This avoids Target plant being same as Source
            ## Calculate d for each Infected plant
            d <- sum(sqrt((Coo[iInfected,1:2] - CooSusceptible)^2))
            ## Calculate the dispersal kernel based on each d and then sum over all kernels
            ## Using Neri et al. 2014 normalized exponential kernel
            m <- m + chooseKernel(distance = d, alpha = alpha, kernelFunc = "normalized exponential")
          }
          ## Evaluate if cryptic period is exceeded
          ## Infectious plants become Diseased when time step t exceeds Inf_times + It
          ## Only update Infectious plants that aren't Diseased yet
          if(t >= (infTimeInfected + It[iInfected]) & Disease_times[iInfected] == Tmax){
            if(Disease_times[iInfected] < t){
              warning("Old disease time being updated")
            }
            Disease_times[iInfected] <- runif(1, min = t-1, max = t)
          }
        } # End iInfected loop
        ###############################################################################################
        #### Epsilon calculations
        ## If vectorOverwintering == TRUE and it's the start of the year, external vector infectivity and density....
        ## ...are dependent on in-field infectivity and density at the end of the previous year
        if(vectorOverwintering == TRUE & y > 1){
          ## For each year after the first ...
          ## ... in-field infectivity at end of previous year contributes to external infectivity but decays exponentially as season progresses
          kappa_e <- kappa_e0 + VOdecay(kappa_b_ym1, t)
        } else {
          kappa_e <- kappa_e0
        }
        ## Immigrating vector density as damped sine wave
        rho_et <- vectorDensity(A = A, lambda_osc = lambda_osc, phi = phi, base = base, time = t)
        ## Save rho_et
        rho_etVec[t] <- rho_et
        ## Calculate epsilonti
        epsilonti <- eta*kappa_e*rho_et*epsilonK[iiSusceptible]
        ## Save epsilonti
        epsilonMatrix[iiSusceptible, t] <- epsilonti
        ############################################################################################################
        #### Beta calculations
        ## if() statement for start of the year
        if(t == 1){
          rho_bt_tm1 <- 0 # Initial in-field density 
          rho_et_tm1 <- 0 # Initial immigrating density 
          kappa_b <- 0 # Initial in-field infectivity
        } else {
          ## Calculating in-field infectivity
          Infecteds_tm1 <- sum(Inf_times < t) # Number of Infectious plants at t-1
          Diseased_tm1 <- sum(Disease_times < t) # Number of Diseased plants at t-1
          kappa_b <- (aI*Infecteds_tm1 + aD*Diseased_tm1)/numPlants # Time-dependent in-field infectivity
          ## Getting vector densities (in-field and immigrating) from previous time step
          rho_bt_tm1 <- rho_btVec[t-1]
          rho_et_tm1 <- rho_etVec[t-1]
        }
        kappaMatrix[t,y] <- kappa_b + kappa_e # Save kappa_b + kappa_e
        ## In-field vector density
        rho_btVec[t] <- (1 - muv)*(rho_bt_tm1 + rho_et_tm1) # Calculate in-field vector density
        beta_t <- eta*kappa_b*rho_btVec[t] # Calculate beta_t
        betaMatrix[iiSusceptible, t] <- beta_t*m # Save beta_ti
        ## Cumulative force of infection for iiSusceptible plant
        lambda[iiSusceptible] <- lambda[iiSusceptible] + beta_t*m + epsilonti
        ## If cumulative lambda (i.e., infection pressure) exceeds Q (i.e., Sellke threshold) plant becomes infected,
        ## then receives an infection time in range [t - 1, t]
        if(lambda[iiSusceptible] >= Q[iiSusceptible]){
          Exp_times[iiSusceptible, y] <- runif(1, min = t-1, max = t)
        }
      } # End iiSusceptible loop
      ########################################################################################
      #### Moving from Exposed to Infectious compartments
      ## Select indices of plants that are in Exposed compartment (but not Infectious compartment yet)
      Exposed_ind <- which(Exp_times[,y] < Tmax & Inf_times == Tmax)
      ## Evaluate if latent period is exceeded
      ## Exposed plants become Infectious when time step t exceeds Exp_times + Et 
      for(iExposed in Exposed_ind){
        if(Inf_times[iExposed] < t){
          warning("Old infection time incorrectly being updated")
        }
        if(t >= Exp_times[iExposed, y] + Et[iExposed]){
          Inf_times[iExposed] <- runif(1, min = t-1, max = t)
        }
      } # End iExposed loop
    } # End t loop
    ###########################################################################################
    #### Winter processes
    ###########################################################################################
    ## If simulation is in it's last year, skip winter processes
    if(y < numYears){
      #### Host recovery
      ## Base host recovery on Exposure time, rather than Infection time
      ## Get all plants that became Exposed; differs from Exposed_ind above
      expIndices <- which(Exp_times[,y] < Tmax)
      recoveryStatus <- probRecovery <- rep(0, length(numPlants))
      for(iExposed in expIndices){
        ## hostRecovery function is based on Feil et al. 2003 equation based on day of year; need to convert to week time step
        probRecovery[iExposed] <- hostRecovery(Exp_times[iExposed, y]*7, b = b)
        recoveryStatus[iExposed] <- rbinom(1, 1, probRecovery[iExposed])
        ## Make sure recoveryStatus is either 0 or 1
        if(recoveryStatus[iExposed] != 0 & recoveryStatus[iExposed] != 1){
          warning("Something is wrong with recoveryStatus")
        }
        ## For chronic infected plants, set Exp_times for next year to 0
        if(recoveryStatus[iExposed] == 0){
          Exp_times[iExposed, y+1] <- 0
        }
      } # End Recovery iExposed loop 
      ## Check that recovery probability calculations are working correctly
      # recoveryCheck <- cbind(Exp_times[expIndices, c(y, y+1)], probRecovery[expIndices], recoveryStatus[expIndices])
      ###########################################################################################
      #### Roguing (removal) and replanting of diseased hosts
      ## Assume that only Diseased (i.e., symptomatic) plants are eligible for removal...
      ## ... and that all Diseased plants have equal probability of removal.
      ## Assume that managers have are limited by the total number of plants they can rogue and replant, rather than a fraction.
      diseaseIndices <- which(Disease_times < Tmax)
      if(length(diseaseIndices) <= numPlantsRogued){
        ## If all diseased plants are to be rogued
        Exp_times[diseaseIndices, y+1] <- Tmax
      } else {
        ## Get random sub-sample of plants to rogue and replace
        roguedPlantIndices <- sample(diseaseIndices, size = numPlantsRogued, replace = FALSE)
        ## Reset Exposed times to Tmax, making them Susceptible again
        Exp_times[roguedPlantIndices, y+1] <- Tmax
      } # End Roguing 
    } # End Winter if() statement
  } # End y loop
  resultsList <- list(Coo = Coo, # Plant coordinates
                      Exp_times = Exp_times, # Vector of Exposure times
                      Inf_times = Inf_times, # Vector of Infection times
                      Disease_times = Disease_times, # Vector of Disease times
                      rho_btVec = rho_btVec, # In-field vector density in last year
                      kappaMatrix = kappaMatrix, # Natural infectivity for every time step and year
                      epsilonMatrix = epsilonMatrix, # Epsilon for every plant and time step in last year
                      betaMatrix = betaMatrix) # Beta for every time step in last year
  return(resultsList)
}



