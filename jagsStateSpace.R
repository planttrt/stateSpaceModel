library(rjags)
##### test with simulations
source('stateSpaceModel.R')

ssSim <- ssSimulations(nSites = 5, nTSet = 10, beta = c(1.,  2),
                       sig = .1, tau = .01, 
                       plotFlag = T, miss = .7, TRUNC = F)


Fore <- which(is.na(ssSim$connect[,1]))
Back <- which(is.na(ssSim$connect[,2]))
Both <- which(!rowSums(is.na(ssSim$connect)))

model <- ifelse(ssSim$TRUNC, 'modelSS.trunc.bugs','modelSS.bugs')

jags <- jags.model(model,
                   data = list('x' = ssSim$x,
                               'z' = ssSim$z,
                               'N' = nrow(ssSim$x),
                               'p'= ncol(ssSim$x),
                               #'Truncated' = ssSim$TRUNC,
                               connectFore=Fore,
                               connectBoth=Both,
                               connectBack=Back),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 2000)

tmp <- jags.samples(jags,c('y','beta', 'sigma', 'sigma2'),2000 )



print(tmp)

par(mfrow=c(3,1))
plot(tmp$beta)
abline(h=ssSim$beta, col='red')

plot(NA, ylim = c(0,1), xlim=c(1, length(tmp$sigma)))
lines(tmp$sigma, type = 'l')
abline(h=sqrt(ssSim$sig), col='red')
lines(tmp$sigma2)
abline(h=sqrt(ssSim$tau), col='red')


beta <- t(apply(tmp$beta, c(1,2), mean))
#hist(beta[-c(1:500),2], breaks = 20)
x <- ssSim$x
y <- ssSim$y
sigma <- t(apply(tmp$sigma, 1:2, mean))
tau <- t(apply(tmp$sigma2, 1:2, mean))
plotType <- 'all'
nPlot =10


latent <- apply(tmp$y, c(1,2), mean)
confInterval <- apply(latent, 1, quantile, probs=c(.025,.5,.975))
cols <- rep('blue', length(ssSim$y))
n <- length(ssSim$y)

cols[ssSim$wNA] <- 'black'

plot(ssSim$y, col=cols)
lines(confInterval[2,])
polygon(c(1:n, n:1), c(confInterval[1,], rev(confInterval[3,])), col='#88888888')
lines(confInterval[1,], col='grey')
lines(confInterval[3,], col='grey')
segments(x0 = 1:n, y0 = confInterval[1,], 
         x1 = 1:n, y1 = confInterval[3,],
         col = 'grey')
points(ssSim$y, col=cols)


