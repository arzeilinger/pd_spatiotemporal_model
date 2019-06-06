#### PROCESSING AND ANALYZING MCMC OUTPUT FROM CLUSTER

rm(list = ls())

my.packages <- c("dplyr", "tidyr", "ggplot2", "coda", "data.table")
lapply(my.packages, require, character.only = TRUE)


source("R_functions/plotting_functions.R")


#### Read in raw MCMC samples from cluster run
samples <- readRDS("output/raw_mcmc_samples.rds")
str(samples)
samples1 <- samples[[1]]


## Getting just the model parameters
paramcols <- which(attr(samples1, "dimnames")[[2]] == "alpha" | attr(samples1, "dimnames")[[2]] == "beta" | attr(samples1, "dimnames")[[2]] == "epsilon")
params1 <- samples1[,paramcols]
## Make paramList a list of MCMC samples of just the three parameters
paramList <- list(params1, samples[[2]][,paramcols])

pdf("output/trace_density_plots_parameters_Adrakey_simulation_cluster.pdf")
  plot(params1)
dev.off()

samplesPlot(samples1, 'epsilon', width = 4, height = 2, traceplot = TRUE, densityplot = FALSE)



## Examining both chains 

chainsPlot(paramList, nrows = 2, buffer.right = 1)

## Gelman-Rubin diagnostic
(grDiagnostics <- gelman.diag(paramList))

## Effective sample size
(essTable <- round(cbind(
  length = apply(params1, 2, length),
  ESS    = effectiveSize(params1)
)))




## Summary of posterior distributions
res <- round(cbind(
  `cilower` = apply(samples1, 2, function(x) quantile(x, 0.025)),
  `mean`    = apply(samples1, 2, mean),
  `median`  = apply(samples1, 2, median),
  `ciupper` = apply(samples1, 2, function(x) quantile(x, 0.975))
), 8)

tail(res)
