R scripts for stochastic spatio-temporal epidemic models

Based on the model of Adrakey et al. 2017 Proc. Royal Society

nimble_model directory contains nimble model and script to simulate infection dynamics. Primary script of interest for running nimble version of model is in nimble_Adrakey_model_perry.R.

R_functions directory contains separate scripts for simulateDiseaseSpread, getTimeIntervals, dispersal kernels, and MCMC samplers (sampler_infection.R).
