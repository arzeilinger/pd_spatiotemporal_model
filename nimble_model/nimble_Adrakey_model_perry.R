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
    #P <- 0
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

## Testing dDiseaseSpread code
# dDiseaseSpread(x = Inf_indices,
#                alpha = 0.08, 
#                beta = 4e-05,
#                epsilon = 6e-04, 
#                Coo = Coo[1:numPlants, 1:2], 
#                Inf_times = Inf_times[1:numPlants], 
#                Tmax = Tmax,
#                numPlants = numPlants, 
#                numInfections = numInfections)

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
Tmax <- 360
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
tail(Y)

#### Set up data set from simulated infection dynamics
simulationResults <- readRDS("output/simulation_results_list_Adrakey_setup.rds")
str(simulationResults)
## Make all elements of the saved list available to global environment
list2env(simulationResults, globalenv())
## Change name of time_intervals
Y <- time_intervals
Coo <- as.matrix(Coo)


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

## MCMC sampler_infection
source("R_functions/sampler_infection.R")

## Version of sampler_infection without the printed text (except when an error is retuned)
source("R_functions/sampler_infection_quiet.R")

## version of buildMCMC for debugging
source("nimble_model/buildMCMC_debug.R")

## Set initial parameter values to "true values" from Adrakey et al. Suppl Mat
Cmodel$alpha <- 0.08
Cmodel$beta <- 7e-06
Cmodel$epsilon <- 5e-05

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
#MCMCconf$removeSamplers(c("alpha", "beta"))
#MCMC <- buildMCMC_debug(MCMCconf)
MCMC <- buildMCMC(MCMCconf)
MCMCconf$printSamplers()
MCMCconf$printMonitors()

## So that set.seed is the same for compiled and uncompiled runs
seed <- 1

## Compile MCMC for running in C++
Cmcmc <- compileNimble(MCMC, project = Rmodel, resetFunctions = TRUE)
 
# ## Run compiled MCMC
# set.seed(seed)
# Cmcmc$run(niter = 100, reset = TRUE)
# 
# 
# ## Run uncompiled MCMC
# set.seed(seed)
# MCMC$run(niter = 100)


####################################################################################################
#### Large MCMC run on Adrakey simulated data
niter <- 20000
nburnin <- 1000
thin <- 1
nchains <- 1
## Check returned no. of samples
floor((niter-nburnin)/thin)

system.time(samples <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, thin = thin, 
                               nchains = nchains, setSeed = seed, samplesAsCodaMCMC = TRUE))
saveRDS(samples, file = "output/raw_mcmc_samples_my_simulated_data_2019-10-17.rds")

## 100000 iterations took 9.6 hours

## Summary of posterior distributions
res <- round(cbind(
  `cilower` = apply(samples, 2, function(x) quantile(x, 0.025)),
  `mean`    = apply(samples, 2, mean),
  `median`  = apply(samples, 2, median),
  `ciupper` = apply(samples, 2, function(x) quantile(x, 0.975))
  ), 8)

tail(res)


#### Examining mixing
library(coda)
source("R_functions/plotting_functions.R")
## Getting just the model parameters
paramcols <- which(attr(samples, "dimnames")[[2]] == "alpha" | attr(samples, "dimnames")[[2]] == "beta" | attr(samples, "dimnames")[[2]] == "epsilon")
params <- samples[,paramcols]

#pdf("output/trace_density_plots_parameters_Adrakey_simulation.pdf")
  plot(as.mcmc(params))
#dev.off()

samplesPlot(samples, c('epsilon'), width = 4, height = 2, traceplot = TRUE, densityplot = FALSE)


#### Looking at infection times
meanInfTimes <- res[grep("Inf_times", row.names(res)), "mean"]
hist(meanInfTimes, breaks = seq(0,Tmax,by=30))
hist(Inf_times, breaks = seq(0,490,by=30))
## Mean epidemic size at 340
sum(meanInfTimes <= 340) # From the MCMC
sum(Inf_times <= 340) # From the data
## Mean epidemic size at Tmax = 460
sum(meanInfTimes < 460)
summary(meanInfTimes)


sum(Inf_times <= 130)
sort(Inf_times[Inf_times <= 130])

sum(Inf_times <= 260)
