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
