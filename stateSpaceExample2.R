

source('stateSpaceModel2.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/dataViz.R')
source('~/Projects/procVisData/colorSet.R')
source('~/Projects/procVisData/colorProc.R')
source('~/Projects/procVisData/bayesianFunctions.R')

set.seed(4)

ssSim <- ssSimulations(nSites = 1, #number of sites
                       nTSet = 20, #number of Time steps
                       beta = c(2, 1, -1), #beta coefficients
                       sig = .02, # process error
                       tau = .2, # observation error
                       plotFlag = T, # whether plot the data or not
                       miss = 0.0, #fraction of missing data
                       TRUNC = F, # trunctaed model or not
                       nonLinear =F, #nonlinear term,
                       ymax=8
)




ssOut <- stateSpace(x = ssSim$x, #predictors
                    z = ssSim$z,#response
                    tauPriorMean = ssSim$tau, 
                    tauPriorStrength = .05, 
                    sigPriorMean = ssSim$sig,
                    sigPriorStrength = 2,
                    connect = ssSim$connect, # connectivity of time data
                    ng = 3000, #number of gibbs sampling
                    burnin = 1000, #numbe of burning
                    MPstep = .1,  # MetroPolis-Hasting initial step size  
                    nAccept = 30,  # number of steps to check the acceptance rate
                    storeLatent = T, # whether store latent state or not
                    TRUNC = ssSim$TRUNC, # whether trucated or not
                    nonLinear = F
)

out <- stateSpaceJags(ssSim$x, ssSim$z, connect = ssSim$connect, nGibbs = 1000)
str(ssOut)

summary(ssOut$bgibbs)
