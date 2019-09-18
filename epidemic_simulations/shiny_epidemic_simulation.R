#### Epidemic simulation for Cal Academy NightLife

rm(list = ls())

my.packages <- c("ggplot2", "dplyr", "tidyr", "raster")
lapply(my.packages, require, character.only = TRUE)

#### Load simulation function and dispersal kernel functions
source("R_functions/simulateDiseaseSpread.R")
source("R_functions/dispersal_kernel_functions.R")
source("R_functions/getTimeIntervals.R")


#### Visualize normalized dispersal kernel
dvec <- c(1:20)
normkout <- normalizedKernel(dvec, alpha = 5)
tiff("output/normalized_exponential_dispersal_kernel_functional_plot.tiff")
  plot(x = dvec, y = normkout, type = "l", lwd = 2, xlab = "Distance", ylab = "K(d, alpha)")
dev.off()

#### Simulation for visualization
#### Initial values setup

# ## Parameter values
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


#### Run simulateDiseaseSpread
testSpread <- simulateDiseaseSpread(alpha, beta, epsilon, Tmax, Coo)
(Inf_times <- testSpread$Inf_times)


###################################################################################################################
#### Plotting infection dynamics
## Produce maps of infection status based on different time points t

## Vector of time points to plot
tcuts <- floor(seq(5,100,length.out=9))


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
plot(diseaseRasterStack, legend = FALSE)


pdf("output/simulated_disease_spread_raster_stack.pdf")
plot(diseaseRasterStack)
dev.off()

#### Make a GIF of the raster stack
giftitle <- rep("Vineyard disease spread", dim(diseaseRasterStack)[3])
raster::animate(diseaseRasterStack, legend = FALSE, main = giftitle, n = 5)

saveGIF(raster::animate(diseaseRasterStack, legend = FALSE, main = giftitle, n = 5), 
        movie.name = "disease_spread.gif")
## Check if it worked correctly
ani.options("convert")
## If it worked correctly, this should return NULL
