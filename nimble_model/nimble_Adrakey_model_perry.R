require(nimble)

#### Simple exponential kernel as nimbleFunction
kernelnf <- nimbleFunction(
    run = function(distance = double(0), alpha = double(0)) {
        return(exp(-(distance / alpha)))
        returnType(double())        
    }
    )

#### Normalized exponential kernel as nimbleFunction
normalizedKernelnf <- nimbleFunction(
    run = function(distance = double(0), alpha = double(0)) {
        return((1/(2*3.141593*distance*alpha))*exp(-(distance / alpha)))
        returnType(double())        
    }
    )


#### nimbleFunction for Adrakey complete data likelihood equation
## Based on Suppl. Material Eqn 1 and likelihood.h C++ code

dDiseaseSpread <- nimbleFunction(
  run = function(x = double(1),
    alpha = double(0), 
    beta = double(0), 
    epsilon = double(0), 
    Coo = double(2), 
    Inf_times = double(1), # We assume these are sorted.
    Tmax = double(0),
    numPlants = double(0),
    numInfections = double(0), 
    log = integer(0, default = 0)) {
    returnType(double(0))
    ## x is Inf_indices
    # numInfections <- length(Inf_indices)
    # numPlants <- nrow(Coo)
    ## Initialize S, for "summation", what's in the big parentheses of equation (1)
    S <- 0
    ## Probabilities of 
    for(iiSourcePlant in 1:numInfections) {
      iSourcePlant <- x[iiSourcePlant]
      infTimeSource <- Inf_times[iSourcePlant]
      CooSource <- Coo[iSourcePlant,1:2]
      for(iTargetPlant in 1:numPlants) {
        #infTimeTarget is either the infection time, or, if never infected, Tmax.
        infTimeTarget <- Inf_times[iTargetPlant]
        if(infTimeTarget > infTimeSource) { ## This also avoids Target plant being same as Source
          d <- sqrt(sum((CooSource - Coo[iTargetPlant,1:2])^2))
          # if(d==0) {
          #   message("case 1")
          #   browser()
          # }
          S <- S + beta * normalizedKernelnf(distance = d, alpha = alpha) * (infTimeTarget - infTimeSource)
        }
      }
    }
    ## We'll move the epsilon term here:
    S <- S + epsilon * sum(Inf_times[1:numPlants])
    ## Need only those plants that were infected by t[obs] and need to sort them by infection time
    ## Need function to sort Inf_times and Inf_indices within nimble
    ## We think that the Adrakey code missed epsilon for the first infected plant
    P <- log(epsilon) ## log force of infection for first infected plant
    for(iTargetPlant in 2:numInfections) {
      m <- 0
      for(iSourcePlant in 1:(iTargetPlant-1)) {
        d <- sqrt(sum((Coo[iSourcePlant,1:2] - Coo[iTargetPlant,1:2])^2)) 
        # if(d==0) {
        #   message("case 2")
        #   browser()
        # }
        m <- m + beta * normalizedKernelnf(distance = d, alpha = alpha)
      }
      P <- P + log(m + epsilon)
    }
    logProb <- P-S
    if(log) return(logProb)
    else return(exp(logProb))
  })


#### nimbleFunction for observations
dDiseaseObs <- nimbleFunction(
  run = function(x = double(0), 
                 Inf_indices = double(1),
                 Inf_times = double(1),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    ## For now, make function return 1 if log = TRUE, 0 if log = FALSE
    if(log) return(0)
    else return(1)
})


#### Set up data set from simulated infection dynamics
simulationResults <- readRDS("output/simulation_results_list.rds")
str(simulationResults)
## Make all elements of the saved list available to global environment
list2env(simulationResults, globalenv())
## Change name of time_intervals
Y <- time_intervals

code <- nimbleCode({
  alpha ~ dunif(0, 1000)
  beta ~ dunif(0, 1000)
  epsilon ~ dunif(0, 1000)
  latent_period <- 100
  ## numInfections and Inf_times are important but will be modified only by custom samplers,
  ## so they do not need their own declarations.
  Inf_indices[1:numPlants] ~ dDiseaseSpread(alpha = alpha, 
                                            beta = beta,
                                            epsilon = epsilon, 
                                            Coo = Coo[1:numPlants, 1:2], 
                                            Inf_times = Inf_times[1:numPlants], 
                                            Tmax = Tmax,
                                            numPlants = numPlants, 
                                            numInfections = numInfections)
  Ydummy ~ dDiseaseObs(Inf_indices[1:numPlants],
                       Inf_times[1:numPlants])
})

constants <- list(## We don't want Tmax in constants because we need it to be a node,
  ## but a scalar constant does not become a node in the model
                  numPlants = numPlants,
                  Coo = Coo)

data <- list(Ydummy = 1)

inits <- list(alpha = 0.1,
              beta = 0.1,
              epsilon = 0.1,
              numInfections = numInfections,
              Tmax = Tmax,
              Inf_times = Inf_times,
              Inf_indices = Inf_indices)

Rmodel <- nimbleModel(code = code, constants = constants, data = data, inits = inits)
Cmodel <- compileNimble(Rmodel, resetFunctions = TRUE)   


sampler_infection <- nimbleFunction(
  name = 'sampler_infection',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## target should be the Inf_indices nodes
    ## time_intervals, latent_period_node, Inf_times_node, and numInfections_node are required
    Inf_indices_node <- target
    time_intervals <- if(!is.null(control$time_intervals))
      control$time_intervals
    else stop("must provide time_intervals to sampler_infection")
    latent_period_node <- if(!is.null(control$latent_period_node))
      control$latent_period_node
    else stop("must provide latent_period_node to sampler_infection")
    Inf_times_node <- if(!is.null(control$Inf_times_node))
      control$Inf_times_node
    else stop("must provide Inf_times_node to sampler_infection")
    numInfections_node <- if(!is.null(control$numInfections_node))
      control$numInfections_node
    else stop("must provide numInfections_node to sampler_infection")
     Tmax_node <- if(!is.null(control$Tmax_node))
      control$Tmax_node
    else stop("must provide Tmax_node to sampler_infection")

     ## optional control list entries
     trials <- if(!is.null(control$trials))  control$trials  else 10
     
     numPlants <- length(model[[Inf_times_node]])
     logProb_Inf_indices_node <- paste0("logProb_", Inf_indices_node)
     # Conventional way to determine what needs to be calculated:   
     calcNodes <- model$getDependencies(target)
     # In this case, we are going to make strong assumptions about model structure anyway,
     # So we could create dependencies by hand.
  },
  run = function() {
    for(k in 1:trials) {
        print("on trial ", k)
        iPlant <- ceiling(runif(1,0,numPlants))
      if(time_intervals[iPlant, 1] != model[[Tmax_node]]) {
        ## known symptomatic
        if(time_intervals[iPlant, 2] > time_intervals[iPlant, 1])
          update_infection_time(iPlant, time_intervals[iPlant, 1], time_intervals[iPlant, 2])
      } else {
        if(model[[Inf_times_node]][iPlant] < model[[Tmax_node]]) {
          ## never observed as symptomatic but has an imputed infection time
          ## within Dt (latent period) of Tmax.
          if(runif(1) < 0.5) {
            ## propose that it stays infected and gets a new time
            update_infection_time(iPlant, model[[Tmax_node]] - model[[latent_period_node]], model[[Tmax_node]])
          } else {
            ## propose that it isn't infected.
            propose_remove_infected(iPlant)
          }
        } else {
          ## never observed as symptomatic and does not have an imputed infection time,
          ## so propose that it does.
          propose_add_infected(iPlant)
        }
      }
    }
  },
 methods = list(
     check_system = function() {
         print("in check_system")
     ok <- TRUE
     current_numInfections <- model[[numInfections_node]]
     if(current_numInfections > 1) {
       for(i in 1:(current_numInfections-1)) {
         if(model[[Inf_times_node]][ model[[Inf_indices_node]][i] ] >=
            model[[Inf_times_node]][ model[[Inf_indices_node]][i+1] ])
           ok <- FALSE
       }
     }
     if(model[[Inf_times_node]][ model[[Inf_indices_node]][current_numInfections] ] == model[[Tmax_node]])
       ok <- FALSE
     if(current_numInfections < numPlants) {
       if(!all(model[[Inf_times_node]][ model[[Inf_indices_node]][(current_numInfections+1):numPlants]] ==
               model[[Tmax_node]]))
         ok <- FALSE
     }
                  print("leaving check_system")
     return(ok)
     returnType(logical())
   },
   propose_add_infected = function(iPlant = double()) {
     ## put iPlant into the numInfections+1 position on Inf_indices
     ## Then do a regular update_infection_time by with different acceptance probability
       ## via the logProb_RJ argument
       print("entering propose_add_infected")
     if(!check_system()) {
       print(model[[Inf_times_node]])
       print(model[[Inf_indices_node]])
       stop("Found the system out of valid state entering propose_add_infected")
       #browser()
     }
     current_numInfections <- model[[numInfections_node]]
     if(current_numInfections < numPlants) {
       proposal_numInfections <- current_numInfections+1
       model[[Inf_indices_node]][proposal_numInfections] <<- iPlant
       model[[numInfections_node]] <<- proposal_numInfections
       startPossibleTime <- model[[Tmax_node]] - model[[latent_period_node]]
       endPossibleTime <- model[[Tmax_node]]
       logProb_RJ_contributions <- log(0.5) + log(endPossibleTime - startPossibleTime)
       jumped <- update_infection_time(iPlant, startPossibleTime, endPossibleTime,
                                       logProb_RJ_contributions)
       ## Inf_indices and Inf_times will have been updated by update_infection_time.
       ## It doesn't matter that iPlant was inserted beyond Inf_indices[numInfections], even if rejected.
       ## But numInfections needs to be updated here
       if(jumped) {
         mvSaved[numInfections_node, 1] <<- proposal_numInfections
       } else {
         model[[numInfections_node]] <<- current_numInfections
       }
     }
     if(!check_system()) {
       print(model[[Inf_times_node]])
       print(model[[Inf_indices_node]])
       stop("Found the system out of valid state exiting propose_add_infected")
       #browser()
     }
       print("leaving propose_add_infected")
   },
   propose_remove_infected = function(iPlant = double()) {
     ## This can work similarly to propose_add_infection, 
     ## but the startPossibleTime and endPossibleTime passed to update_infection_time will both be Tmax
       ## so the proposal is forced to be Tmax
              print("entering propose_remove_infected")
     if(!check_system()) {
       print(model[[Inf_times_node]])
       print(model[[Inf_indices_node]])
       stop("Found the system out of valid state entering propose_remove_infected")
       #browser()
     }

     current_numInfections <- model[[numInfections_node]]
     if(current_numInfections > 0) {
       # These start and end possible times are for the logProb_RJ_contributions
       startPossibleTime <- model[[Tmax_node]] - model[[latent_period_node]]
       endPossibleTime <- model[[Tmax_node]]
       logProb_RJ_contributions <- - log(0.5) - log(endPossibleTime - startPossibleTime)
       ## We provide endPossibleTime for both start and end, so that the proposal will definitely be at the end.
       jumped <- update_infection_time(iPlant, endPossibleTime, endPossibleTime,
                                       logProb_RJ_contributions)
       ## Inf_indices and Inf_times will have been updated by update_infection_time.
       ## It doesn't matter that iPlant was inserted beyond Inf_indices[numInfections], even if rejected.
       ## But numInfections needs to be updated here
       if(jumped) {
         mvSaved[numInfections_node, 1] <<- current_numInfections - 1
       }
       if(!check_system()) {
         print(model[[Inf_times_node]])
         print(model[[Inf_indices_node]])
         stop("Found the system out of valid state exiting propose_remove_infected")
         #browser()
       }

     }
       print("exiting propose_remove_infected")
   },
   update_infection_time = function(iPlant = double(),
                                    startPossibleTime = double(),
                                    endPossibleTime = double(),
                                    logProb_RJ_contributions = double(0, default = 0)) {
       returnType(logical(0))
                            print("entering update_infection_time")
     if(!check_system()) {
       print(model[[Inf_times_node]])
       print(model[[Inf_indices_node]])
       stop("Found the system out of valid state entering update_infection_time")
       #browser()
     }
    currentLogProb <- model$getLogProb(calcNodes)
    current_iSorted <- 1
    numInfections <- model[[numInfections_node]]
    while(model[[Inf_indices_node]][current_iSorted] != iPlant
          & current_iSorted <= numInfections) {
      current_iSorted <- current_iSorted + 1
    }
    if(current_iSorted > numInfections)
      stop("Problem in update_known_symptomatic: current_iSorted > numInfections")
    
    current <- model[[Inf_times_node]][iPlant]
    ## startPossibleTime might equal endPossibleTime
    proposal <- startPossibleTime + runif(1) * (endPossibleTime - startPossibleTime);

    if(proposal < current) {
      if(current_iSorted == 1) {
        proposal_iSorted <- 1
      } else {
        proposal_iSorted <- current_iSorted
        done <- FALSE
        while(proposal_iSorted > 1 & !done)
          if(model[[Inf_times_node]][model[[Inf_indices_node]][proposal_iSorted - 1] ] > proposal) {
            proposal_iSorted <- proposal_iSorted - 1
          } else {
            done <- TRUE
          }
      }
    } else {
      proposal_iSorted <- current_iSorted
      done <- FALSE
      while(proposal_iSorted < numInfections & ! done) {
        if(model[[Inf_times_node]][model[[Inf_indices_node]][proposal_iSorted + 1] ] < proposal) {
          proposal_iSorted <- proposal_iSorted + 1
        } else {
          done <- TRUE
        }
      }
    }
    if(proposal_iSorted != current_iSorted) {
      ## If the plant is moved in the indices -- based on the proposal -- this routine shifts the indices of other plants accordingly
      if(proposal < current) {
        moveStartFrom <- proposal_iSorted
        moveEndFrom <- current_iSorted-1
        moveStartTo <- moveStartFrom + 1
        moveEndTo <- moveEndFrom + 1
      } else {
        moveStartFrom <- current_iSorted + 1
        moveEndFrom <- proposal_iSorted
        moveStartTo <- moveStartFrom - 1
        moveEndTo <- moveEndFrom - 1
      }
      model[[Inf_indices_node]][moveStartTo:moveEndTo] <<- model[[Inf_indices_node]][moveStartFrom:moveEndFrom]
    }
    model[[Inf_indices_node]][proposal_iSorted] <<- iPlant
    model[[Inf_times_node]][iPlant] <<- proposal
    if(!check_system()) {
      print(model[[Inf_times_node]])
      print(model[[Inf_indices_node]])
      stop("Found the system out of valid state in update_infection_time before calculating proposal logProb.")
      #browser()
    }
    proposalLogProb <- model$calculate(calcNodes)
    logAcceptProb <- proposalLogProb - currentLogProb + logProb_RJ_contributions
    outcome <- decide(logAcceptProb)
                                        #if(decide(logAcceptProb)) {
       print("outcome = ", outcome)
       print(moveStartTo, " ", moveEndTo, " ", moveStartFrom, " ", moveEndFrom, " ", proposal_iSorted, current_iSorted)
    if(outcome) {
      ## accept
      ## The updates after accepting or rejecting are customized to this particular model.
      ## We need to be sure that if we change the model, we check if this code needs modification.
      ## This code should move all modified parts of the model into the correspoding "1" row of mvSaved.
      if(proposal_iSorted != current_iSorted) {
        mvSaved[Inf_indices_node, 1][moveStartTo:moveEndTo] <<- model[[Inf_indices_node]][moveStartTo:moveEndTo] 
      }
      mvSaved[Inf_times_node, 1][iPlant] <<- proposal
      mvSaved[Inf_indices_node, 1][proposal_iSorted] <<- iPlant
      mvSaved[logProb_Inf_indices_node, 1][1] <<- model[[logProb_Inf_indices_node]][1]
      jumped <- TRUE
    } else {
      ## reject
      ## This code should restore all modified parts of the model to what they were upon entry to this method.
        if(proposal_iSorted != current_iSorted) {
            model[[Inf_indices_node]][moveStartTo:moveEndTo] <<- mvSaved[Inf_indices_node, 1][moveStartTo:moveEndTo]
##        model[[Inf_indices_node]][moveStartFrom:moveEndFrom] <<- model[[Inf_indices_node]][moveStartTo:moveEndTo]
      }
      model[[Inf_times_node]][iPlant] <<- current
      model[[Inf_indices_node]][proposal_iSorted] <<- mvSaved[Inf_indices_node, 1][proposal_iSorted]
      model[[logProb_Inf_indices_node]][1] <<- mvSaved[logProb_Inf_indices_node, 1][1]
      jumped <- FALSE
    }
    if(!check_system()) {
       print(model[[Inf_times_node]])
       print(model[[Inf_indices_node]])
       stop("Found the system out of valid state exiting update_infection_time.")
      #browser()
    }
       print("exiting update_infection_time")
    return(jumped)
   },
  reset = function() {
  })
)

MCMCconf <- configureMCMC(Rmodel, nodes = c("alpha", "beta", "epsilon"))
MCMCconf$addSampler(type = 'sampler_infection', 
                    target = 'Inf_indices',
                    control = list(
                      time_intervals = Y,
                      Inf_times_node = 'Inf_times',
                      Inf_indices_node = 'Inf_indices',
                      numInfections_node = 'numInfections',
                      latent_period_node = 'latent_period',
                      Tmax_node = 'Tmax'
                    ))
MCMCconf$addMonitors('Inf_times')
source("nimble_model/buildMCMC_debug.R")
MCMC <- buildMCMC_debug(MCMCconf)
MCMCconf$printSamplers()


set.seed(123)

## Compile MCMC for running in C++
Cmcmc <- compileNimble(MCMC, project = Rmodel, resetFunctions = TRUE)

Cmcmc$run(niter = 100)

## Run compiled MCMC
samples <- runMCMC(Cmcmc, niter = 100, nchains = 1)

## Run uncompiled MCMC
samples <- runMCMC(MCMC, niter = 100, nchains = 1)

