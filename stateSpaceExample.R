

source('stateSpaceModel.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/dataViz.R')
source('~/Projects/procVisData/bayesianFunctions.R')
par(mfrow=c(2,2))

# set.seed(2016)

ssSim <- ssSimulations(nSites = 8, #number of sites
                       nTSet = 10, #number of Time steps
                       beta = c(3,2, -2, -1), #beta coefficients
                       sig = .1, # process error
                       tau = .3, # observation error
                       plotFlag = F, # whether plot the data or not
                       miss = .3, #fraction of missing data
                       TRUNC = F # trunctaed model or not
)

ssOut <- stateSpace(x = ssSim$x, #predictors
                    z = ssSim$z,#response
                    tauPriorMean = NULL, 
                    tauPriorStrength = NULL, 
                    sigPriorMean = ssSim$sig,
                    sigPriorStrength = 2,
                    connect = ssSim$connect, # connectivity of time data
                    ng = 10000, #number of gibbs sampling
                    burnin = 5000, #numbe of burning
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
#ssGibbs <- ssGibbs[-(1:ssOut$nBurnin),]
trueValues <- c(ssSim$beta, ssSim$sig, ssSim$tau)
yRange <- range(ssGibbs[-(1:ssOut$nBurnin),])

# plotGibbsDensity(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, 
#           trueValues = trueValues,
#           yRange = c(0,15))

plotGibbsChains(outGibbs =  ssGibbs, colPalette = colorRampPalette(topo.colors(ncol(ssGibbs))),# burnin = ssOut$nBurnin, 
                 trueValues = trueValues)

plotGibbsBoxplots(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, sort = F, textAdj = 0, sigPlot = F,
                  trueValues = trueValues, vert = T, rangeY = range(c(trueValues, yRange)))

stateSpaceTemporalPost(ssSim$x, ssSim$y, 
                       beta = ssOut$bgibbs,
                       sigma = ssOut$sgibbs, 
                       startPoints = ssSim$startPoints, 
                       nTrends = 100)

points(ssSim$wNA, ssSim$y[ssSim$wNA], col='blue')

stateSpacePlot(y=ssSim$y, 
               yGibbs=  ssOut$latentGibbs, 
               wMissing= ssSim$wNA )
