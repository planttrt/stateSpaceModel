source('stateSpaceAuxFunctions.R')
source('stateSpaceModel.R')

ssSim <- ssSimulations(nSites = 3, #number of sites
                       nTSet = 20, #number of Time steps
                       p = 3, # number of predictors
                       sig = .1, # process error
                       tau = .01, # observation error
                       plotFlag = T, # whether plot the data or not
                       TRUNC = T # trunctaed model or not 
)

ssOut <- stateSpace(x = ssSim$x, #predictors
                    z = ssSim$z,#response
                    tauPriorMean = ssSim$tau, 
                    tauPriorStrength = 100, 
                    sigPriorMean = NULL,
                    sigPriorStrength = NULL,
                    connect = ssSim$connect, # connectivity of time data
                    ng = 2000, #number of gibbs sampling
                    burnin = 500, #numbe of burning
                    MPstep = .1,  # MetroPolis-Hasting initial step size  
                    nAccept = 30,  # number of steps to check the acceptance rate
                    TRUNC = ssSim$TRUNC # whether trucated or not
)


plot(rowMeans(ssOut$mpSteps), type = 'l')
plot(rowMeans(ssOut$acceptGibbs), type = 'l')

plot(ssOut$bgibbs[,1]); abline(h=ssSim$beta[1], col='red')
plot(ssOut$bgibbs[,2]); abline(h=ssSim$beta[2], col='red')
plot(ssOut$bgibbs[,3]); abline(h=ssSim$beta[3], col='red')

plot(ssOut$sgibbs); abline(h=ssSim$sig, col='red')
plot(ssOut$tgibbs); abline(h=ssSim$tau, col='red')

summary(ssOut$bgibbs[-(1:ssOut$nBurnin),])
ssSim$beta

plot(ssSim$y, ssOut$latentMean, 
     xlab = 'Latent state - observed', 
     ylab = 'Latent state - predicted')

abline(0,1, col='#3366EE')

plot(ssSim$z, colMeans(ssOut$zPredGibbs[-c(1:ssOut$nBurnin),]), 
     xlab = 'Response - observed', 
     ylab = 'Response - predicted')

abline(0,1, col='#BB3333')

