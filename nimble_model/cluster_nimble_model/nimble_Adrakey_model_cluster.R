#### Running Adrakey spatiotemporal epidemic model on the cluster

## Need to specify library location for nimble because of issues installing it on the cluster
library(foreach)
library(doParallel)
library(nimble, lib.loc = '~/R/tmp')
library(coda, lib.loc = '~/R/tmp')

## Specify SLURM environmental variables for parallelization
nCores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
registerDoParallel(nCores)


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



#####################################################################################################
#### Data sets


#### For sampler_infection to run correctly data sets need to be structured correctly, and slightly differently from Adrakey et al.:
## Tmax represents the last observation time (syn. with t[obs] in Adrakey et al.)
## time_interval needs to be a 2-column matrix, 
## with the first column as the last time the plant was observed asymptomatic
## and the second column as the first time the plant was observed diseased.
## If a plant was never observed diseased, both columns of time_interval coded as Tmax and Inf_times == Tmax

#### This is the data set supplied by Adrakey
adcoo <- read.table("Adrakey2017_code/Coo.txt")
str(adcoo)
adIndx_ind <- read.table("Adrakey2017_code/Indx_ind.txt")
str(adIndx_ind)
adInf_tim <- read.table("Adrakey2017_code/Inf_tim.txt")
str(adInf_tim)
adtime_interval <- read.table("Adrakey2017_code/time_interval.txt")
str(adtime_interval)

Coo <- as.matrix(adcoo)
Inf_indices <- as.numeric(adIndx_ind[,1])
Tmax <- 460
## Correct Tmax entries in Inf_times, now that Tmax == 460 (or Adrakey's t[obs])
Inf_times <- adInf_tim[,1]
Inf_times[Inf_times > Tmax] <- Tmax
numPlants <- nrow(adcoo)
numInfections <- sum(Inf_times < Tmax)
## Correct and rename adtime_interval
Y <- as.matrix(adtime_interval)
for(i in 1:nrow(Y)){
  if(Y[i,1] == 0 & Y[i,2] == 0){
    Y[i,] <- c(Tmax, Tmax)
  }
}

#### Set up data set from simulated infection dynamics
# simulationResults <- readRDS("output/simulation_results_list.rds")
# str(simulationResults)
# ## Make all elements of the saved list available to global environment
# list2env(simulationResults, globalenv())
# ## Change name of time_intervals
# Y <- time_intervals


####################################################################################
#### Nimble model code

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

## Set parameters to "true values" from Adrakey et al. Suppl Mat
Cmodel$alpha <- 0.08
Cmodel$beta <- 7e-06
Cmodel$epsilon <- 5e-05

## MCMC sampler_infection
#source("R_functions/sampler_infection.R")

## Version of sampler_infection without the printed text (except when an error is retuned)
source("R_functions/sampler_infection_quiet.R")

## version of buildMCMC for debugging
#source("nimble_model/buildMCMC_debug.R")

MCMCconf <- configureMCMC(Rmodel, nodes = c("alpha", "beta", "epsilon"))
MCMCconf$addSampler(type = 'sampler_infection_quiet', 
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
#MCMC <- buildMCMC_debug(MCMCconf)
MCMC <- buildMCMC(MCMCconf)
MCMCconf$printSamplers()

## So that set.seed is the same for compiled and uncompiled runs
#seed <- 1

## Compile MCMC for running in C++
Cmcmc <- compileNimble(MCMC, MCMC2, project = Rmodel, resetFunctions = TRUE)

## Run compiled MCMC
#set.seed(seed)
#Cmcmc$MCMC$run(niter = 100, reset = TRUE)


####################################################################################################
#### Large MCMC run on Adrakey simulated data
niter <- 40000
nburnin <- 1000
thin <- 1
nchains <- 2
## Check returned no. of samples
#print(floor((niter-nburnin)/thin))

ti <- Sys.time()

samples <- foreach(mcmcseed = 1:2) %dopar% {
	runMCMC(Cmcmc, niter = niter, nburnin = nburnin, thin = thin,
		nchains = nchains, setSeed = mcmcseed, samplesAsCodaMCMC = TRUE)
}

tf <- Sys.time()
print(tf - ti)

saveRDS(samples, "output/raw_mcmc_samples_test_2019-09.rds")

samples1 <- samples[[1]]

## Summary of posterior distributions
res <- round(cbind(
  `cilower` = apply(samples1, 2, function(x) quantile(x, 0.025)),
  `mean`    = apply(samples1, 2, mean),
  `median`  = apply(samples1, 2, median),
  `ciupper` = apply(samples1, 2, function(x) quantile(x, 0.975))
  ), 8)

print(tail(res))


#### Examining mixing
## Getting just the model parameters
#paramcols <- which(attr(samples1, "dimnames")[[2]] == "alpha" | attr(samples1, "dimnames")[[2]] == "beta" | attr(samples1, "dimnames")[[2]] == "epsilon")
# params <- samples1[,paramcols]
# 
# pdf("output/trace_density_plots_parameters_Adrakey_simulation.pdf")
#   plot(as.mcmc(params))
# dev.off()

print("SUCCESS!!")
