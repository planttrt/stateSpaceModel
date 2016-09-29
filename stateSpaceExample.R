

source('stateSpaceModel.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/dataViz.R')
source('~/Projects/procVisData/bayesianFunctions.R')
par(mfrow=c(3,2))

set.seed(2016)

ssSim <- ssSimulations(nSites = 4, #number of sites
                       nTSet = 10, #number of Time steps
                       beta = c(1., -.3, -.4, .2), #beta coefficients
                       sig = .2, # process error
                       tau = .01, # observation error
                       plotFlag = T, # whether plot the data or not
                       miss = .3, #fraction of missing data
                       TRUNC = F # trunctaed model or not
)

ssOut <- stateSpace(x = ssSim$x, #predictors
                    z = ssSim$z,#response
                    tauPriorMean = ssSim$tau, 
                    tauPriorStrength = 1, 
                    sigPriorMean = NULL,
                    sigPriorStrength = NULL,
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

plotGibbsDensity(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, 
          trueValues = c(ssSim$beta, ssSim$sig, ssSim$tau),
          yRange = c(0,15))

plotGibbsChains(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, 
                 trueValues = c(ssSim$beta, ssSim$sig, ssSim$tau))

plotGibbsBoxplots(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, sort = F)

stateSpaceTemporalPost(ssSim$x, ssSim$y, 
                       beta = ssOut$bgibbs,
                       sigma = ssOut$sgibbs, 
                       startPoints = ssSim$startPoints, 
                       nTrends = 100)

points(ssSim$wNA, ssSim$y[ssSim$wNA], col='blue')

stateSpacePlot(y=ssSim$y, 
               yGibbs=  ssOut$latentGibbs, 
               wMissing= ssSim$wNA )
