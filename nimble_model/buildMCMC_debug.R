## copied from nimble source code and modified with output messages
buildMCMC_debug <- nimbleFunction(
    name = 'MCMC_debug',
    setup = function(conf, ...) {
    	if(inherits(conf, 'modelBaseClass'))   conf <- configureMCMC(conf, ...)
    	else if(!inherits(conf, 'MCMCconf')) stop('conf must either be a nimbleModel or a MCMCconf object (created by configureMCMC(...) )')
        dotdotdotArgs <- list(...)
        enableWAICargument <- if(!is.null(dotdotdotArgs$enableWAIC)) dotdotdotArgs$enableWAIC else nimbleOptions('enableWAIC')    ## accept enableWAIC argument regardless
        model <- conf$model
        my_initializeModel <- initializeModel(model)
        mvSaved <- modelValues(model)
        samplerFunctions <- nimbleFunctionList(sampler_BASE)
        for(i in seq_along(conf$samplerConfs))
            samplerFunctions[[i]] <- conf$samplerConfs[[i]]$buildSampler(model=model, mvSaved=mvSaved)
        samplerExecutionOrderFromConfPlusTwoZeros <- c(conf$samplerExecutionOrder, 0, 0)  ## establish as a vector
        monitors  <- mcmc_processMonitorNames(model, conf$monitors)
        monitors2 <- mcmc_processMonitorNames(model, conf$monitors2)
        thinFromConfVec <- c(conf$thin, conf$thin2)  ## vector
        thinToUseVec <- c(0, 0)                      ## vector, needs to member data
        mvSamplesConf  <- conf$getMvSamplesConf(1)
        mvSamples2Conf <- conf$getMvSamplesConf(2)
        mvSamples <- modelValues(mvSamplesConf)
        mvSamples2 <- modelValues(mvSamples2Conf)
        samplerTimes <- c(0,0) ## establish as a vector
        progressBarLength <- 52  ## multiples of 4 only
        progressBarDefaultSetting <- getNimbleOption('MCMCprogressBar')
        ## WAIC setup:
        dataNodes <- model$getNodeNames(dataOnly = TRUE)
        dataNodeLength <- length(dataNodes)
        sampledNodes <- model$getVarNames(includeLogProb = FALSE, nodes = monitors)
        sampledNodes <- sampledNodes[sampledNodes %in% model$getVarNames(includeLogProb = FALSE)]
        paramDeps <- model$getDependencies(sampledNodes, self = FALSE, downstream = TRUE)
        allVarsIncludingLogProbs <- model$getVarNames(includeLogProb = TRUE)
        enableWAIC <- enableWAICargument || conf$enableWAIC   ## enableWAIC comes from MCMC configuration, or from argument to buildMCMC
        if(enableWAIC) {
            if(dataNodeLength == 0)   stop('WAIC cannot be calculated, as no data nodes were detected in the model.')
            mcmc_checkWAICmonitors(model = model, monitors = sampledNodes, dataNodes = dataNodes)
        }
    },

    run = function(
        niter                 = integer(),
        reset                 = logical(default = TRUE),
        time                  = logical(default = FALSE),
        progressBar           = logical(default = TRUE),
        ## reinstate samplerExecutionOrder as a runtime argument, once we support non-scalar default values for runtime arguments:
        ##samplerExecutionOrder = integer(1, default = -1)
        nburnin               = integer(default =  0),
        thin                  = integer(default = -1),
        thin2                 = integer(default = -1)) {
        if(nburnin < 0)   stop('cannot specify nburnin < 0')
        thinToUseVec <<- thinFromConfVec
        if(thin  != -1)   thinToUseVec[1] <<- thin
        if(thin2 != -1)   thinToUseVec[2] <<- thin2
        print("initializing model")
        my_initializeModel$run()
        print("done initializing model. copying to mvSaved")
        nimCopy(from = model, to = mvSaved, row = 1, logProb = TRUE)
                print("done copying to mvSaved")
        if(reset) {
            for(i in seq_along(samplerFunctions))   samplerFunctions[[i]]$reset()
            samplerTimes <<- numeric(length(samplerFunctions) + 1)       ## default inititialization to zero
            mvSamples_copyRow  <- 0
            mvSamples2_copyRow <- 0
        } else {
            if(nburnin != 0)   stop('cannot specify nburnin when using reset = FALSE.')
            if(dim(samplerTimes)[1] != length(samplerFunctions) + 1)   samplerTimes <<- numeric(length(samplerFunctions) + 1)   ## first run: default inititialization to zero
            mvSamples_copyRow  <- getsize(mvSamples)
            mvSamples2_copyRow <- getsize(mvSamples2)
        }
        print("done resetting")
        resize(mvSamples,  mvSamples_copyRow  + floor((niter-nburnin) / thinToUseVec[1]))
        resize(mvSamples2, mvSamples2_copyRow + floor((niter-nburnin) / thinToUseVec[2]))
        ## reinstate samplerExecutionOrder as a runtime argument, once we support non-scalar default values for runtime arguments:
        ##if(dim(samplerExecutionOrder)[1] > 0 & samplerExecutionOrder[1] == -1) {   ## runtime argument samplerExecutionOrder was not provided
        ##    lengthSamplerExecutionOrderFromConf <- dim(samplerExecutionOrderFromConfPlusTwoZeros)[1] - 2
        ##    if(lengthSamplerExecutionOrderFromConf == 0) samplerExecutionOrderToUse <- numeric(0) else samplerExecutionOrderToUse <- samplerExecutionOrderFromConfPlusTwoZeros[1:lengthSamplerExecutionOrderFromConf]
        ##} else {   ## runtime argument samplerExecutionOrder was provided
        ##    samplerExecutionOrderToUse <- samplerExecutionOrder
        ##}
        lengthSamplerExecutionOrderFromConf <- dim(samplerExecutionOrderFromConfPlusTwoZeros)[1] - 2
        if(lengthSamplerExecutionOrderFromConf == 0) samplerExecutionOrderToUse <- numeric(0) else samplerExecutionOrderToUse <- samplerExecutionOrderFromConfPlusTwoZeros[1:lengthSamplerExecutionOrderFromConf]
        if(niter < progressBarLength+3 | !progressBarDefaultSetting) progressBar <- progressBar & 0  ## cheap way to avoid compiler warning
        if(progressBar) { for(iPB1 in 1:4) { cat('|'); for(iPB2 in 1:(progressBarLength/4)) cat('-') }; print('|'); cat('|') }
        progressBarIncrement <- niter/(progressBarLength+3)
        progressBarNext <- progressBarIncrement
        progressBarNextFloor <- floor(progressBarNext)
        if(niter < 1) return()
        print("starting samplers")
        for(iter in 1:niter) {
            checkInterrupt()
            if(time) {
                for(i in seq_along(samplerExecutionOrderToUse)) {
                    print("running sampler", i)
                    ind <- samplerExecutionOrderToUse[i]
                    samplerTimes[ind] <<- samplerTimes[ind] + run.time(samplerFunctions[[ind]]$run())
                    print("done with sampler", i)
                }
            } else {
                for(i in seq_along(samplerExecutionOrderToUse)) {
                    print("running sampler", i)
                    ind <- samplerExecutionOrderToUse[i]
                    samplerFunctions[[ind]]$run()
                    print("done with sampler", i)
                }
            }
            if(iter > nburnin) {
                print("recording params")
                sampleNumber <- iter - nburnin
                if(sampleNumber %% thinToUseVec[1] == 0) {
                    mvSamples_copyRow  <- mvSamples_copyRow  + 1
                    nimCopy(from = model, to = mvSamples,  row = mvSamples_copyRow,  nodes = monitors)
                }
                if(sampleNumber %% thinToUseVec[2] == 0) {
                    mvSamples2_copyRow <- mvSamples2_copyRow + 1
                    nimCopy(from = model, to = mvSamples2, row = mvSamples2_copyRow, nodes = monitors2)
                }
                print("done recording params")
            }
            if(progressBar & (iter == progressBarNextFloor)) {
                cat('-')
                progressBarNext <- progressBarNext + progressBarIncrement
                progressBarNextFloor <- floor(progressBarNext)
            }
        }
        if(progressBar) print('|')
        returnType(void())
    },
    methods = list(
        getTimes = function() {
            returnType(double(1))
            return(samplerTimes[1:(length(samplerTimes)-1)])
        },
        calculateWAIC = function(nburnin = integer(default = 0),
            burnIn = integer(default = 0)) {
            if(!enableWAIC) {
                print('Error: must set enableWAIC = TRUE in buildMCMC. See help(buildMCMC) for additional information.')
                return(NaN)
            }
            if(burnIn != 0) {
                print('Warning: \'burnIn\' argument is deprecated and will not be supported in future versions of NIMBLE. Please use the \'nburnin\' argument instead.')
                ## If nburnin has not been changed, replace with burnIn value
                if(nburnin == 0)   nburnin <- burnIn
            }
            nburninPostThinning <- ceiling(nburnin/thinToUseVec[1])
            numMCMCSamples <- getsize(mvSamples) - nburninPostThinning
            if((numMCMCSamples) < 2) {
                print('Error: need more than one post burn-in MCMC samples')
                return(-Inf)
            }
            logPredProbs <- matrix(nrow = numMCMCSamples, ncol = dataNodeLength)
            logAvgProb <- 0
            pWAIC <- 0
            currentVals <- values(model, allVarsIncludingLogProbs)
            
            for(i in 1:numMCMCSamples) {
                copy(mvSamples, model, nodesTo = sampledNodes, row = i + nburninPostThinning)
                model$simulate(paramDeps)
                model$calculate(dataNodes)
                for(j in 1:dataNodeLength)
                    logPredProbs[i,j] <- model$getLogProb(dataNodes[j])
            }
            for(j in 1:dataNodeLength) {
                maxLogPred <- max(logPredProbs[,j])
                thisDataLogAvgProb <- maxLogPred + log(mean(exp(logPredProbs[,j] - maxLogPred)))
                logAvgProb <- logAvgProb + thisDataLogAvgProb
                pointLogPredVar <- var(logPredProbs[,j])
                pWAIC <- pWAIC + pointLogPredVar
            }
            WAIC <- -2*(logAvgProb - pWAIC)
            values(model, allVarsIncludingLogProbs) <<- currentVals
            if(is.nan(WAIC)) print('WAIC was calculated as NaN.  You may need to add monitors to model latent states, in order for a valid WAIC calculation.')
            returnType(double())
            return(WAIC)
        }),
    where = getLoadingNamespace()
)
