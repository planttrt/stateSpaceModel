

source('stateSpaceModel.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/dataViz.R')
source('~/Projects/procVisData/colorSet.R')
source('~/Projects/procVisData/colorProc.R')
source('~/Projects/procVisData/bayesianFunctions.R')

 set.seed(4)

ssSim <- ssSimulations(nSites = 3, #number of sites
                       nTSet = 12, #number of Time steps
                       beta = c(2, -1, 1), #beta coefficients
                       sig = .02, # process error
                       tau = .2, # observation error
                       plotFlag = F, # whether plot the data or not
                       miss = 0.8, #fraction of missing data
                       TRUNC = F # trunctaed model or not
)

ssOut <- stateSpace(x = ssSim$x, #predictors
                    z = ssSim$z,#response
                    tauPriorMean = ssSim$tau, 
                    tauPriorStrength = .05, 
                    sigPriorMean = ssSim$sig,
                    sigPriorStrength = 2,
                    connect = ssSim$connect, # connectivity of time data
                    ng = 2000, #number of gibbs sampling
                    burnin = 1000, #numbe of burning
                    MPstep = .1,  # MetroPolis-Hasting initial step size  
                    nAccept = 30,  # number of steps to check the acceptance rate
                    storeLatent = T, # whether store latent state or not
                    TRUNC = ssSim$TRUNC # whether trucated or not
)


# plot(rowMeans(ssOut$mpSteps), type = 'l')
# plot(rowMeans(ssOut$acceptGibbs), type = 'l')


ssGibbs <- data.frame(beta=ssOut$bgibbs,
                      sigma=ssOut$sgibbs,
                      tau=ssOut$tgibbs)
ssGibbs.betas <- ssGibbs[,1:(ncol(ssGibbs)-2)]
ssGibbs.errors <- ssGibbs[,(ncol(ssGibbs)-1):ncol(ssGibbs)]

#ssGibbs <- ssGibbs[-(1:ssOut$nBurnin),]
trueValues <- c(ssSim$beta, ssSim$sig, ssSim$tau)
trueValues.betas <- ssSim$beta
trueValues.errors <- c(ssSim$sig, ssSim$tau)
yRange <- range(ssGibbs[-(1:ssOut$nBurnin),])

# plotGibbsDensity(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, 
#           trueValues = trueValues,
#           yRange = c(0,15))


png(filename = paste0('example/post-',ssSim$miss,'.png'), width = 8, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))

plotGibbsChains(outGibbs =  ssGibbs.betas, trueValues = trueValues.betas, 
                title.text = '(a) betas', ylab = '')
plotGibbsChains(outGibbs =  ssGibbs.errors, trueValues = trueValues.errors,
                title.text = '(b) tau and sigma', ylab = '')

plotGibbsBoxplots(outGibbs =  ssGibbs.betas, burnin = ssOut$nBurnin, sort = F, textAdj = 0, sigPlot = F,
                  trueValues = trueValues.betas, vert = T)
mtext('(c) betas', font = 2)

plotGibbsBoxplots(outGibbs =  ssGibbs.errors, burnin = ssOut$nBurnin, sort = F, textAdj = 0, sigPlot = F,
                  trueValues = trueValues.errors, vert = T, rangeY = c(-10*ssSim$sig,5*ssSim$tau))
mtext('(d) tau and sigma', font = 2)

dev.off()


png(filename = paste0('example/pred-',ssSim$miss,'.png'), width = 6, height = 6, units = 'in', res = 300)
par(mfrow=c(2,1), mar=c(4,4,2,1))
stateSpacePlot(y=ssSim$y, 
               startPoints = ssSim$startPoints, 
               yGibbs=  ssOut$latentGibbs, 
               wMissing= ssSim$wNA , yLim = range(ssSim$z, na.rm = T)*1.2)
mtext('(e) Predictions of one step ahead', font = 2)
stateSpaceTemporalPost(ssSim$x, ssSim$z, 
                       beta = ssOut$bgibbs,
                       sigma = ssOut$sgibbs, 
                       startPoints = ssSim$startPoints, 
                       nTrends = 500, 
                       yLim = range(ssSim$z, na.rm = T)*1.2)
mtext('(f) Predictions of the whole time series from the intial point', font = 2)

dev.off()


