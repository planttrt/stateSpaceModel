library(rjags)

source('stateSpaceModel.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/gibbsSamplingPlot.R')

par(mfrow=c(2,2))

set.seed(27708)
ssSim <- ssSimulations(nSites = 1, #number of sites
                       nTSet = 50, #number of Time steps
                       beta = c(.1, -.07), #beta coefficients
                       sig = .1, # process error
                       tau = .01, # observation error
                       plotFlag = T, # whether plot the data or not
                       miss = .1, #fraction of missing data
                       TRUNC = T # trunctaed model or not
)

ssJAGS <- stateSpaceJags(x = ssSim$x, z = ssSim$z, connect = ssSim$connect, 
               nGibbs = 2000, nPost = 2000, TRUNC = ssSim$TRUNC, quiet = T)

plotGibbs(outGibbs =  ssJAGS$chains, burnin = 1000,
          trueValues = c(ssSim$beta, ssSim$sig, ssSim$tau), 
          #yRange = c(0,15),
          plots = 'post')

stateSpaceTemporalPost(ssSim$x, ssSim$y, 
                       beta = ssJAGS$chains[,1:ncol(ssSim$x)],
                       sigma = ssJAGS$chains[,ncol(ssSim$x)+1], 
                       startPoints = ssSim$startPoints, 
                       nTrends = 100)

stateSpacePlot(y=ssSim$y, 
               yGibbs=  t(apply(ssJAGS$rawsamples$y, c(1,2), mean)), 
               wMissing= ssSim$wNA )
