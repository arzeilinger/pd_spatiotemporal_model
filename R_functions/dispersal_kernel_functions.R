#### Dispersal kernel functions for Adrakey spatiotemporal epidemic model

#### Simple exponential kernel
kernel <- function(distance, alpha) {
  exp(-(distance / alpha))
}

#### Normalized exponential kernel
## Derived in Neri et al. 2014
## Used in Adrakey et al. 2017
normalizedKernel <- function(distance, alpha) {
  (1/(2*pi*distance*alpha))*exp(-(distance / alpha))
}

#### Cauchy kernel function
## From Neri et al. 2014 but not normalized
cauchyKernel <- function(distance, alpha){
  2/(pi*alpha*(1+((distance^2)/(alpha^2))))
}

#### Wrapper function to choose among different kernel functions
chooseKernel <- function(distance, alpha, kernelFunc = "exponential"){
  if(kernelFunc == "exponential"){
    kda <- kernel(distance = distance, alpha = alpha)
  } else {
    if(kernelFunc == "normalized exponential"){
      kda <- normalizedKernel(distance = distance, alpha = alpha)
    } else {
      if(kernelFunc == "Cauchy"){
        kda <- cauchyKernel(distance = distance, alpha = alpha)
      }
    }
  }
  return(kda)
}

#### Function to make Coordinates on a grid
makeCoordinates <- function(nrc){
  ## Number of rows and columns are the same -- will always make the plant population a square
  ## Deviating from a square causes a problem with making time slice rasters -- not sure why
  grid <- 1:nrc
  Coo <- matrix(c(rep(grid, nrc), rep(grid, each = nrc)), ncol = 2)
  ## Re-order in descending order of Coo[,2]
  ## Note, the ordering of the coordinates only matters for plotting results as rasters
  ## Raster values are set starting at the top-left corner and going across
  Coo <- Coo[order(Coo[,2], decreasing = TRUE),]
  return(Coo)
}