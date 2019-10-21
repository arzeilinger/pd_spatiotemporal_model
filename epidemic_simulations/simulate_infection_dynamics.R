#### FUNCTION TO SIMULATE EPIDEMIC PROCESS
## Result is a simulated data set from infection dynamices
## According to Adrakey's force of infection equation

rm(list = ls())

my.packages <- c("ggplot2", "dplyr", "tidyr", "raster", "animation")
lapply(my.packages, require, character.only = TRUE)

#### Load simulation function and dispersal kernel functions
source("R_functions/simulateDiseaseSpread.R")
source("R_functions/getTimeIntervals.R")

#### Simulation for visualization
#### Initial values setup

#### Parameter values
alpha <- 0.8 # Dispersal parameter
beta <- 0.3 # Contact rate or secondary spread rate
epsilon <- 0.001 # Primary spread rate

## Number of time steps
Tmax <- 200

#### Plant spatial coordinates
## Number of rows and columns are the same -- will always make the plant population a square
## Deviating from a square causes a problem with making time slice rasters -- not sure why
nrc <- 20
grid <- 1:nrc
Coo <- matrix(c(rep(grid, nrc), rep(grid, each = nrc)), ncol = 2)


#### Run simulateDiseaseSpread
testSpread <- simulateDiseaseSpread(alpha, beta, epsilon, Tmax, Coo)
(Inf_times <- testSpread$Inf_times)
Inf_indices <- testSpread$Inf_indices

#saveRDS(testSpread, "output/test_stochastic_simulation_infection_dynamics.rds")


###################################################################################################################
#### Plotting infection dynamics
## Produce maps of infection status based on different time points t

## Vector of time points to plot
tcuts <- floor(seq(1,Tmax, by = 5))


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


# pdf("output/simulated_disease_spread_raster_stack.pdf")
#   plot(diseaseRasterStack)
# dev.off()


#### Make a GIF of the raster stack
giftitle <- rep("Vineyard disease spread: \ngrey = healthy, green = infected", dim(diseaseRasterStack)[3])
raster::animate(diseaseRasterStack, legend = FALSE, main = giftitle, n = 5)

saveGIF(raster::animate(diseaseRasterStack, legend = FALSE, main = giftitle, n = 5), 
        movie.name = "disease_spread_stochastic_test.gif")
## Check if it worked correctly
ani.options("convert")
## If it worked correctly, this should return NULL


##########################################################################################################
#### Produce "snapshots" of disease observations
## tcuts (observation times) should probably not include Tmax (i.e., the last observation should be < Tmax)
## tcuts needs to start with 0
tslices <- seq(0, Tmax, by = 10)
time_intervals <- getTimeIntervals(Inf_times = Inf_times, tcuts = tslices) ## Doesn't work for some reason, maybe related to running (and stopping) the GIF beforehand?
cbind(time_intervals, Inf_times)

#### Compile all the simulation components
simulationResults <- list(Tmax = Tmax,
                          alpha = alpha,
                          beta = beta,
                          epsilon = epsilon,
                          latent_period = 3,
                          Coo = Coo,
                          numPlants = nrow(Coo),
                          numInfections = sum(Inf_times < Tmax),
                          Inf_times = Inf_times,
                          Inf_indices = Inf_indices)
                          #time_intervals = time_intervals)
saveRDS(simulationResults, "output/simulation_results_list.rds")



####################################################################################################
#### Replicating simulation from Adrakey et al. 2017
#### Warning: simulation takes about 1 hour on Berkeley computing cluster

## Parameter values used in simulations
## Units are included as comments
alphaAd <- 0.08 # km^2
betaAd <- 7e-6 # day^-1 km^2
epsilonAd <- 5e-5 # day^-1

TmaxAd <- 360

## Simulated host population
## N = 1000 plants, uniformly distributed over 0.75 x 0.75 km^2
## Data from Adrakey et al. 2017 supplementary files
adcoo <- read.table("Adrakey2017_code/Coo.txt")
summary(adcoo)


#### Simulate spread from Adrakey parameters
## This takes > 5 hours on my computer; use the cluster
spreadAd <- simulateDiseaseSpread(alpha = alphaAd, beta = betaAd, epsilon = epsilonAd, Tmax = TmaxAd, Coo = adcoo)
head(spreadAd$Inf_times[order(spreadAd$Inf_times)])

saveRDS(spreadAd, "output/simulation_output_adrakey_setup.rds")

## Read in Inf_times and Inf_indices from simulation run on cluster
spreadAd <- readRDS("output/simulation_output_adrakey_setup.rds")
Inf_times <- spreadAd$Inf_times
Inf_indices <- spreadAd$Inf_indices



#### Make time intervals from my simulations
tcuts <- seq(0, 360, by = 30)

time_intervals <- getTimeIntervals(Inf_times = Inf_times, tcuts = tcuts, Tmax = TmaxAd)

checktime_interval <- cbind(time_intervals, Inf_times)
checktime_interval[checktime_interval[,3] < 360,]
checktime_interval

#### Save simulation initial parameters and outputs
simulationResults <- list(Tmax = TmaxAd,
                          alpha = alphaAd,
                          beta = betaAd,
                          epsilon = epsilonAd,
                          Coo = adcoo,
                          numPlants = nrow(adcoo),
                          numInfections = sum(Inf_times < TmaxAd),
                          Inf_times = Inf_times,
                          Inf_indices = Inf_indices,
                          time_intervals = time_intervals)
saveRDS(simulationResults, "output/simulation_results_list_Adrakey_setup.rds")


######################################################################################################
#### Compare my simulation to Infection times provided by Adrakey paper
adInf_tim <- read.table("Adrakey2017_code/Inf_tim.txt")
str(adInf_tim)

## Adrakey's Inf_times
infectedsAd <- adInf_tim$V1[adInf_tim$V1 < 360]
infectionStatusAd <- ifelse(adInf_tim$V1 < 360, 1, 0)
## Symptomatic plants
length(infectedsAd)
hist(infectedsAd)



#### Compare distribution of infection times at t = 130
t1 <- Inf_times[Inf_times <= 130]
t1[order(t1)]
t1[order(t1)] %>% length()

adInf_tim %>% dplyr::filter(V1 <= 130) %>% arrange(V1) %>% t()


#### Total infections
infecteds <- Inf_times[Inf_times < 460]
infectionStatus <- ifelse(Inf_times < 460, 1, 0)
length(infecteds)
hist(infecteds)
# Symptomatic plants
sum(infecteds < 360)

infData <- data.frame("xcoo" = adcoo$V1,
                      "ycoo" = adcoo$V2,
                      "Inf_times" = Inf_times)

#### Make some plots to compare with Adrakey Figure S2
## Infections at t = 130
timeslice <- 130
simulateAd_t130 <- ggplot(data = infData, aes(x=xcoo, y=ycoo)) +
  geom_point(data = infData[infData$Inf_times > timeslice,], colour = "grey", size = 3) +
  geom_point(data = infData[infData$Inf_times <= timeslice-100,], colour = "red", size = 3) +
  geom_point(data = infData[infData$Inf_times <= timeslice & infData$Inf_times > timeslice-100,], colour = "blue", size = 3) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
simulateAd_t130  

ggsave("output/simulation_adrakey_setup_figure_t130.jpg", plot = simulateAd_t130,
       width = 7, height = 7, units = "in")

## Infections at t = 360
timeslice <- 360
simulateAd_final <- ggplot(data = infData, aes(x=xcoo, y=ycoo)) +
  geom_point(colour = "grey", size = 3) +
  geom_point(data = infData[infData$Inf_times <= timeslice-100,], colour = "red", size = 3) +
  geom_point(data = infData[infData$Inf_times < timeslice & infData$Inf_times > timeslice-100,], colour = "blue", size = 3) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank()) 
simulateAd_final  

ggsave("output/simulation_adrakey_setup_figure_time_final.jpg", plot = simulateAd_final,
       width = 7, height = 7, units = "in")




#### Look at time intervals from Adrakey
adtime_interval <- read.table("Adrakey2017_code/time_interval.txt")
str(adtime_interval)
adtobs <- unique(c(adtime_interval$V1, adtime_interval$V2))
adtobs[order(adtobs)]
## Observations were at 30 day intervals up to 360
## Observations stopped at 360 because Adrakey's simulations assumed a latent period of 100 days
adtime_interval2 <- cbind(adtime_interval, adInf_tim)
adtime_interval2[adtime_interval2[,3] < 460, ]


