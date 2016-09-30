

source('stateSpaceModel.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/dataViz.R')
source('~/Projects/procVisData/bayesianFunctions.R')
par(mfrow=c(3,2))

# set.seed(2016)

ssSim <- ssSimulations(nSites = 100, #number of sites
                       nTSet = 20, #number of Time steps
                       beta = c(.7, -.3, .1, -.1, -.4, .2), #beta coefficients
                       sig = .05, # process error
                       tau = .2, # observation error
                       plotFlag = T, # whether plot the data or not
                       miss = .9, #fraction of missing data
                       TRUNC = F # trunctaed model or not
)

ssOut <- stateSpace(x = ssSim$x, #predictors
                    z = ssSim$z,#response
                    tauPriorMean = NULL, 
                    tauPriorStrength = NULL, 
                    sigPriorMean = ssSim$sig,
                    sigPriorStrength = 1,
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
#ssGibbs <- ssGibbs[-(1:ssOut$nBurnin),]

plotGibbsDensity(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, 
          trueValues = c(ssSim$beta, ssSim$sig, ssSim$tau),
          yRange = c(0,15))

plotGibbsChains(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, 
                 trueValues = c(ssSim$beta, ssSim$sig, ssSim$tau))

plotGibbsBoxplots(outGibbs =  ssGibbs, burnin = ssOut$nBurnin, sort = F, textAdj = 0, trueValues = c(ssSim$beta, ssSim$sig, ssSim$tau), vert = F)

stateSpaceTemporalPost(ssSim$x, ssSim$y, 
                       beta = ssOut$bgibbs,
                       sigma = ssOut$sgibbs, 
                       startPoints = ssSim$startPoints, 
                       nTrends = 100)

points(ssSim$wNA, ssSim$y[ssSim$wNA], col='blue')

stateSpacePlot(y=ssSim$y, 
               yGibbs=  ssOut$latentGibbs, 
               wMissing= ssSim$wNA )
