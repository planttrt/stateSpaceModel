##### test with simulations
source('stateSpaceModel.R')

ssSim <- ssSimulations(nSites = 3, nTSet = 50, beta = c(.3, -.2,1),
                       sig = .01, tau = .001, 
                       plotFlag = T, miss = .0, TRUNC = F)


Fore <- which(is.na(ssSim$connect[,1]))
Back <- which(is.na(ssSim$connect[,2]))
Both <- which(!rowSums(is.na(ssSim$connect)))

jags <- jags.model('modelSS.bugs',
                   data = list('x' = ssSim$x,
                               'z' = ssSim$z,
                               'N' = nrow(ssSim$x),
                               'p'= ncol(ssSim$x),
                               connectFore=Fore,
                               connectBoth=Both,
                               connectBack=Back),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 1000)

tmp <- jags.samples(jags,c('beta', 'sigma', 'sigma2'),2000)



print(tmp)

par(mfrow=c(3,1))
plot(tmp$beta)  
abline(h=ssSim$beta, col='red')

plot(tmp$sigma, type = 'l')
abline(h=sqrt(ssSim$sig), col='red')

plot(tmp$sigma2)
abline(h=sqrt(ssSim$tau), col='red')
