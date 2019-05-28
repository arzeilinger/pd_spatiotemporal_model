#### MCMC algorithm sampler_infection

#### nimbleRcall function to insert into sampler_infection

browser2 <- function(){
    browser(text = "opening browser in R")
    NULL
}

check_systemR <- nimbleRcall(prototype = function(){}, 
                             Rfun = 'browser2', returnType = void())


#### MCMC algorithm sampler_infection
#### Based on Adrakey et al. 2017

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
      #print("on trial ", k)
      iPlant <- ceiling(runif(1,0,numPlants))
      #print("iPlant = ", iPlant)
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
      #print("in check_system")
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
      if(!ok)
      check_systemR()
      #print("leaving check_system")
      return(ok)
      returnType(logical())
    },
    propose_add_infected = function(iPlant = double()) {
      ## put iPlant into the numInfections+1 position on Inf_indices
      ## Then do a regular update_infection_time by with different acceptance probability
      ## via the logProb_RJ argument
      #print("entering propose_add_infected")
      if(!check_system()) {
        #browser()
        stop("Found the system out of valid state entering propose_add_infected")
      }
      current_numInfections <- model[[numInfections_node]]
      if(current_numInfections < numPlants) {
        proposal_numInfections <- current_numInfections+1
        replaced_index <- model[[Inf_indices_node]][proposal_numInfections] ## AZ added this 2019-05-23; the plant index that will be replaced by iPlant
        model[[Inf_indices_node]][proposal_numInfections] <<- iPlant
        ## Find the old iPlant in Inf_indices_node and replace it with replaced_index
        ## AZ added this 2019-05-23
        for(i in (proposal_numInfections+1):numPlants){ 
          if(model[[Inf_indices_node]][i] == iPlant) {
            model[[Inf_indices_node]][i] <<- replaced_index
          }
        }
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
          mvSaved[numInfections_node, 1][1] <<- proposal_numInfections ## AZ added [1] 2019-05-29
        } else {
          model[[numInfections_node]] <<- current_numInfections
        }
      }
      if(!check_system()) {
        ## Using nimbleRcall to open browser()
        check_systemR()
        #browser()
        stop("Found the system out of valid state exiting propose_add_infected")
      }
      #print("leaving propose_add_infected")
    },
    propose_remove_infected = function(iPlant = double()) {
      ## This can work similarly to propose_add_infection, 
      ## but the startPossibleTime and endPossibleTime passed to update_infection_time will both be Tmax
      ## so the proposal is forced to be Tmax
      #print("entering propose_remove_infected")
      if(!check_system()) {
        ## Using nimbleRcall to open browser()
        check_systemR()
        #browser()
        stop("Found the system out of valid state entering propose_remove_infected")
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
          model[[numInfections_node]] <<- current_numInfections - 1 ## AZ changed this 2019-05-23
          mvSaved[numInfections_node, 1][1] <<- current_numInfections - 1 ## AZ added [1] 2019-05-29
        }
        if(!check_system()) {
          ## Using nimbleRcall to open browser()
          check_systemR()
          #browser()
          stop("Found the system out of valid state exiting propose_remove_infected")
        }
        
      }
      #print("exiting propose_remove_infected")
    },
    update_infection_time = function(iPlant = double(),
                                     startPossibleTime = double(),
                                     endPossibleTime = double(),
                                     logProb_RJ_contributions = double(0, default = 0)) {
      returnType(logical(0))
      #print("entering update_infection_time")
      crazyCase <- iPlant == 254 | iPlant == 87
      #if(crazyCase) { 
        #print("In crazy case with iPlant == ", iPlant)
        #print("startPossibleTime == ", startPossibleTime)
        #print("endPossibleTime == ", endPossibleTime)
        #print("model[[numInfections_node]] == ", model[[numInfections_node]])
        #print("Inf_times of infected == ", model[[Inf_times_node]][ model[[Inf_indices_node]][ 1:model[[numInfections_node]] ] ])
      #}
      # if(!check_system()) {
      #   ## Using nimbleRcall to open #browser()
      #   check_systemR()
      #   stop("Found the system out of valid state entering update_infection_time")
      # }
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

      #if(crazyCase) print("current_iSorted = ", current_iSorted, " current = ", current, " proposal = ", proposal)
      
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
      # if(crazyCase) {
      #    print("proposal_iSorted = ", proposal_iSorted)
      #    #check_systemR()
      # }
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
        temp <- model[[Inf_indices_node]][moveStartFrom:moveEndFrom]
        model[[Inf_indices_node]][moveStartTo:moveEndTo] <<- temp
      }
      model[[Inf_indices_node]][proposal_iSorted] <<- iPlant
      model[[Inf_times_node]][iPlant] <<- proposal
      ## AZ changed this 2019-05-23
      ## propose_remove_infected necessarily hits this check_system because model[[numInfections_node]] hasn't been updated yet
      ## this extra if() statement goes around this check_system if we're in propose_remove_infected
      if(startPossibleTime != endPossibleTime){  
        if(!check_system()) {
          print("Proposed new index for iPlant == ", model[[Inf_indices_node]][proposal_iSorted])
          print("Proposed new infection time for iPlant == ", model[[Inf_times_node]][iPlant])
          ## Using nimbleRcall to open #browser()
          check_systemR()
          #browser()
          stop("Found the system out of valid state in update_infection_time before calculating proposal logProb.")
        } 
      }
      proposalLogProb <- model$calculate(calcNodes)
      logAcceptProb <- proposalLogProb - currentLogProb + logProb_RJ_contributions
      outcome <- decide(logAcceptProb)
      #if(decide(logAcceptProb)) {
      #print("outcome = ", outcome)
      #print(moveStartTo, " ", moveEndTo, " ", moveStartFrom, " ", moveEndFrom, " ", proposal_iSorted, current_iSorted)
      if(outcome) {
        ## accept
        ## The updates after accepting or rejecting are customized to this particular model.
        ## We need to be sure that if we change the model, we check if this code needs modification.
        ## This code should move all modified parts of the model into the correspoding "1" row of mvSaved.
          if(proposal_iSorted != current_iSorted) {
            temp <- model[[Inf_indices_node]][moveStartTo:moveEndTo] 
            mvSaved[Inf_indices_node, 1][moveStartTo:moveEndTo] <<- temp 
        }
        mvSaved[Inf_times_node, 1][iPlant] <<- proposal
        #mvSaved[Inf_indices_node, 1][proposal_iSorted] <<- iPlant
        mvSaved[Inf_indices_node, 1] <<- model[[Inf_indices_node]] ## AZ changed this 2019-05-23
        mvSaved[logProb_Inf_indices_node, 1][1] <<- model[[logProb_Inf_indices_node]][1]
        jumped <- TRUE
      } else {
        ## reject
        ## This code should restore all modified parts of the model to what they were upon entry to this method.
          if(proposal_iSorted != current_iSorted) {
            temp <-  mvSaved[Inf_indices_node, 1][moveStartTo:moveEndTo]
            model[[Inf_indices_node]][moveStartTo:moveEndTo] <<- temp
          ##        model[[Inf_indices_node]][moveStartFrom:moveEndFrom] <<- model[[Inf_indices_node]][moveStartTo:moveEndTo]
          }
        model[[numInfections_node]] <<- mvSaved[numInfections_node, 1][1] ## PdV added [1], but I don't see why this line is needed.  I don't see this having been modified previously.  It looks like the add-sampler and remove-sampler steps manage model[[numInfections_node]] when necessary.
        model[[Inf_times_node]][iPlant] <<- current
        model[[Inf_indices_node]] <<- mvSaved[Inf_indices_node, 1] ## AZ changed this 2019-05-23 ## This modification will copy the entire vector instead of the single value. Is it necessary?
        #model[[Inf_indices_node]][proposal_iSorted] <<- mvSaved[Inf_indices_node, 1][proposal_iSorted]
        model[[logProb_Inf_indices_node]][1] <<- mvSaved[logProb_Inf_indices_node, 1][1]
        jumped <- FALSE
      }
      ## AZ changed this 2019-05-23
      ## propose_remove_infected necessarily hits this check_system because model[[numInfections_node]] hasn't been updated yet
      ## this extra if() statement goes around this check_system if we're in propose_remove_infected
      # if(crazyCase) { 
      #   print("In crazy case with iPlant == ", iPlant)
      #   print("startPossibleTime == ", startPossibleTime)
      #   print("endPossibleTime == ", endPossibleTime)
      #   print("model[[numInfections_node]] == ", model[[numInfections_node]])
      #   print("Inf_times of infected == ", model[[Inf_times_node]][ model[[Inf_indices_node]][ 1:model[[numInfections_node]] ] ])
      #   print("jumped == ", jumped)
      #   #stop("Stopped at crazyCase")
      # }
      if(startPossibleTime != endPossibleTime){  
        if(!check_system()) {
          ## Using nimbleRcall to open #browser()
          check_systemR()
          #browser()
          stop("Found the system out of valid state exiting update_infection_time.")
        }
      }
      #print("exiting update_infection_time")
      return(jumped)
    },
    reset = function() {
    })
)
