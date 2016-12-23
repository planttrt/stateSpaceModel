source('stateSpaceModel.R')
source('stateSpacePlot.R')
source('~/Projects/procVisData/dataViz.R')
source('~/Projects/procVisData/colorSet.R')
source('~/Projects/procVisData/colorProc.R')
source('~/Projects/procVisData/bayesianFunctions.R')

miss <-0
for(miss in c(0, .1, .3, .5))
{
  
  set.seed(4)
  ssSim <- ssSimulationsMax(nSites = 4, #number of sites
                            nTSet = 20, #number of Time steps
                            beta = c(1, 2), #beta coefficients
                            sig = .03, # process error
                            tau = .3, # observation error
                            plotFlag = T, # whether plot the data or not
                            miss = miss, #fraction of missing data
                            ymax = c(9,5,7, 3) #maximum of saturation trajectory
  )
  
  ssOut <- stateSpace.Max(x = ssSim$x, #predictors  
                          nGibbs = 20000,
                          nBurnin = 10000,
                          z = ssSim$z,#response
                          connect = ssSim$connect, # connectivity of time data
                          quiet=T)
  
  summ <- getGibbsSummaryMaxModel(ssOut, burnin = 1)
  
  png(paste0('example/example.',sprintf(fmt = '%1.1f',ssSim$miss),'.png'), width = 6, height = 8, units = 'in', res = 300)
  
  par(mar=c(2,1,1,0), oma=c(3,3,1,2), cex.axis=1.5)
  layout(matrix(c(1:6,rep(7,4)), nrow = 5, ncol = 2, byrow = T))
  
  ssSim$errors <- cbind(ssSim$sig, ssSim$tau)
  ssSim$betas <- ssSim$beta
  
  for(var in c('betas','ymax','errors')){
    
    ssGibbs <- summ[[which(names(summ)==var)]]
    trueValues <- ssSim[[which(names(ssSim)==var)]]
    
    yRange <- range(trueValues, ssGibbs)
    
    tmp <- boxplot(ssGibbs, names = NA, ylim=yRange, outline = F)
    text(x = 1:ncol(ssGibbs), y = yRange[1]- .15*diff(range(yRange)),
         labels = colnames(ssGibbs), srt=0, xpd=T, col='#666666', cex=1.5)
    segments(x0 = (1:ncol(ssGibbs))-.4, x1 = (1:ncol(ssGibbs))+.4, y0 = trueValues, y1=trueValues, col='red', lwd=2)
    
    par(yaxt='n')
    plotGibbsChains(outGibbs =  ssGibbs, yRange = yRange, txtSize = 1.5,
                    trueValues = trueValues, txtCol = 'black',txtAdj = 0.5,
                    title.text = '', ylab = '', xlab = '')
    par(yaxt='s')
  }
  
  
  yLim <- range(ssSim$z, na.rm = T)*1.2
  stateSpaceTemporalPost(ssSim$x, ssSim$y, xlab = '',ylab = '',
                         nonLinear = T,
                         ymax = summ$ymax,
                         beta = summ$betas,
                         sigma = summ$errors$sigma, 
                         tau = summ$errors$tau,
                         plotZ = T,
                         startPoints = ssSim$startPoints, 
                         nTrends = 500, 
                         yLim = yLim, connectDots = F )
  
  ssSimPlot(z = ssSim$z, connect = ssSim$connect, ylim = yLim, 
            add = T, col = 'blue', pch = 2, lwd = 1.5)
  
  mtext(text = 'time steps', font = 2, side = 1, cex = 1.5, line = 3)
  mtext(text = 'y', font = 2, side = 2, cex = 1.5, line = 2)
  legend('topright', legend = c('predictions','mean','95%-ile','latent states','observations'),
         lty=c(1,1,2,1,1), col=c('grey','black','black','black','blue'),
         pch=c(NA, NA, NA, 1,2), bty='n', cex=1.3, lwd=2)
  dev.off()
  
}