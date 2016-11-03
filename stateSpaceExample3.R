source('stateSpaceModel.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/dataViz.R')
source('~/Projects/procVisData/colorSet.R')
source('~/Projects/procVisData/colorProc.R')
source('~/Projects/procVisData/bayesianFunctions.R')
miss <-0
for(miss in c(0,.1, .3, .5, .7))
  {
  
set.seed(4)
ssSim <- ssSimulationsMax(nSites = 4, #number of sites
                          nTSet = 20, #number of Time steps
                          beta = c(2, 0, 1), #beta coefficients
                          sig = .1, # process error
                          tau = .4, # observation error
                          plotFlag = T, # whether plot the data or not
                          miss = miss, #fraction of missing data
                          ymax = c(10,5,9,2) #maximum of saturation trajectory
                          )

ssOut <- stateSpace.Max(x = ssSim$x, #predictors  
                          nGibbs = 10000,
                            z = ssSim$z,#response
                            connect = ssSim$connect, # connectivity of time data
                            quiet=T)

ng <- nrow(ssOut$chains)
bunrnt <- (1:floor(ng/2))


ssGibbs.betas <- ssOut$chains[,1:ncol(ssSim$x)]
ssGibbs.ymax <- ssOut$chains[,(ncol(ssSim$x)+1):(ncol(ssOut$chains)-2)]
ssGibbs.errors <- ssOut$chains[,(ncol(ssOut$chains)-1):ncol(ssOut$chains)]



png(paste0('example/post.',ssSim$miss,'.png'), width = 6, height = 6, units = 'in', res = 300)
par(mfrow=c(3,2), mar=c(2,3,1,2))

## betas

ssGibbs <- ssGibbs.betas
trueValues <- ssSim$beta
yRange <- range(trueValues, ssGibbs[-bunrnt,])
tmp <- boxplot(ssGibbs[-bunrnt,], names = NA, ylim=yRange)
text(x = 1:ncol(ssGibbs), y = yRange[1]- .15*diff(range(yRange)), labels = colnames(ssGibbs), srt=0, xpd=T, col='#666666')
segments(x0 = (1:ncol(ssGibbs))-.4, x1 = (1:ncol(ssGibbs))+.4, y0 = trueValues, y1=trueValues, col='red', lwd=2)

plotGibbsChains(outGibbs =  ssGibbs, yRange = yRange,
                trueValues = trueValues, 
                title.text = '', ylab = '', xlab = '')

## ymax

ssGibbs <- ssGibbs.ymax
trueValues <- ssSim$ymax
yRange <- range(trueValues, ssGibbs[-bunrnt,])
tmp <- boxplot(ssGibbs[-bunrnt,], names = NA, ylim=yRange)
text(x = 1:ncol(ssGibbs), y = yRange[1]- .15*diff(range(yRange)), labels = colnames(ssGibbs), srt=0, xpd=T, col='#666666')
segments(x0 = (1:ncol(ssGibbs))-.4, x1 = (1:ncol(ssGibbs))+.4, y0 = trueValues, y1=trueValues, col='red', lwd=2)

plotGibbsChains(outGibbs =  ssGibbs, yRange = yRange,
                trueValues = trueValues, 
                title.text = '', ylab = '', xlab = '')

## errors

ssGibbs <- ssGibbs.errors
trueValues <- c(ssSim$sig, ssSim$tau) 
yRange <- range(trueValues, ssGibbs[-bunrnt,])
tmp <- boxplot(ssGibbs[-bunrnt,], names = NA, ylim=yRange)
text(x = 1:ncol(ssGibbs), y = yRange[1]- .15*diff(range(yRange)), labels = colnames(ssGibbs), srt=0, xpd=T, col='#666666')
segments(x0 = (1:ncol(ssGibbs))-.4, x1 = (1:ncol(ssGibbs))+.4, y0 = trueValues, y1=trueValues, col='red', lwd=2)

plotGibbsChains(outGibbs =  ssGibbs, yRange = yRange,
                trueValues = trueValues, 
                title.text = '', ylab = '', xlab = '')


dev.off()

png(paste0('example/pred.',ssSim$miss,'.png'), width = 10, height = 5, units = 'in', res = 300)
par(mar=c(4,4,1,1))
yLim <- range(ssSim$z, na.rm = T)*1.2
stateSpaceTemporalPost(ssSim$x, ssSim$y, nonLinear = T,
                       ymax = as.matrix(ssGibbs.ymax),
                       beta = ssGibbs.betas,
                       sigma = ssGibbs.errors$sigma, 
                       tau = ssGibbs.errors$tau, plotZ = T,
                       startPoints = ssSim$startPoints, 
                       nTrends = 500, 
                       yLim = yLim, connectDots = T )

ssSimPlot(z = ssSim$z, connect = ssSim$connect, ylim = yLim, 
          add = T, col = 'blue', pch = 2, lwd = 1.5)

legend('topright', legend = c('predictions','mean','95%-ile','latent states','observations'),
       lty=c(1,1,2,1,1), col=c('grey','black','black','black','blue'),
       pch=c(NA, NA, NA, 1,2), bty='n', cex=1.3, lwd=2)
dev.off()

}