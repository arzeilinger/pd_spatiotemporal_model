#### FUNCTION TO SIMULATE EPIDEMIC PROCESS
## Result is a simulated data set from infection dynamices
## According to Adrakey's force of infection equation

rm(list = ls())

my.packages <- c("ggplot2", "dplyr", "tidyr", "raster")
lapply(my.packages, require, character.only = TRUE)

#### Kernel functions
## Simple exponential kernel
kernel <- function(distance, alpha) {
  exp(-(distance / alpha))
}

## Normalized exponential kernel
## Derived in Neri et al. 2014
## Used in Adrakey et al. 2017
normalizedKernel <- function(distance, alpha) {
  (1/(2*pi*distance*alpha))*exp(-(distance / alpha))
}


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
  initial_infected_ind <- sample(1:nrow(Coo), 1)
  ## Update Inf_times and Inf_indiceswith initial infected
  Inf_times[initial_infected_ind] <- runif(1)
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


#### Initial values setup

## Parameter values used in Adrakey et al. 2017 simulations
## Units are included as comments
alpha <- 0.08 # km^2
beta <- 7e-6 # day^-1 km^2
epsilon <- 5e-5 # day^-1

## Alternative parameter values
alpha <- 0.08
beta <- 0.5
epsilon <- 0.01

## Other values
Tmax <- 100
## Number of rows and columns -- will always make the plant population a square
## Deviating from a square causes a problem with making time slice rasters -- not sure why
nrc <- 20
grid <- 1:nrc
Coo <- matrix(c(rep(grid, nrc), rep(grid, each = nrc)), ncol = 2)


## Run simulateDiseaseSpread
testSpread <- simulateDiseaseSpread(alpha, beta, epsilon, Tmax, Coo)
(Inf_times <- testSpread$Inf_times)



###################################################################################################################
#### Plotting infection dynamics
## Produce maps of infection status based on different time points t

## Vector of time points to plot
tcuts <- floor(seq(1,100,length.out=9))


#### Produce raster stack of infection status at a series of time points
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

## Empty raster stack
diseaseRasterList <- vector("list", length(tcuts))

for(i in 1:length(tcuts)){
  ## Produce vector of infection status: 1 = infected, 0 = susceptible at a given time point t
  infectionStatus.i <- ifelse(Inf_times < tcuts[i], 1, 0)
  raster.i <- raster(nrows = CooNrows, ymn = CooYmn, ymx = CooYmx,
                     ncols = CooNcols, xmn = CooXmn, xmx = CooXmx,
                     vals = infectionStatus.i)
  diseaseRasterList[[i]] <- raster.i
}

diseaseRasterStack <- stack(diseaseRasterList)
plot(diseaseRasterStack)




pdf("output/simulated_disease_spread_raster_stack.pdf")                                          
  plot(diseaseRasterStack)
dev.off()




# infectionStatus <- ifelse(Inf_times < tcuts[2], 1, 0)
# 
# ## Make Coo a data.frame for plotting
# diseaseData <- data.frame("cooRows" = Coo[,1],
#                           "cooColumns" = Coo[,2],
#                           "infection" = as.integer(infectionStatus))
# 
# diseaseRaster <- with(diseaseData, raster(nrows = length(unique(cooRows)), 
#                                           ncols = length(unique(cooColumns)),
#                                           ymn = min(cooRows), ymx = max(cooRows),
#                                           xmn = min(cooColumns), xmx = max(cooColumns),
#                                           vals = infection))
# pdf("output/simulated_disease_spread_raster.pdf")                                          
#   plot(diseaseRaster)
# dev.off()

# diseasePlot <- ggplot(data = diseaseData, aes(x = columns, y = rows)) +
#   geom_point(aes(colour = infection), size = 6) +
#   theme_bw() + 
#   theme(axis.line = element_line(colour = "black"),
#         axis.text = element_text(size=14),
#         axis.title = element_text(size=16),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black"),
#         panel.background = element_blank()) 
# 
# diseasePlot
# 
