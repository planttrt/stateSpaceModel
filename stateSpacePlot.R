source('~/Projects/procVisData/auxFunctions.R')
stateSpaceTemporalPost <- function(x, y, beta, t=1:length(y),
                                   xLim=NULL, yLim=NULL, 
                                   sigma=0,
                                   plotFlag =T,
                                   nTrends=10,
                                   startPoints=1,
                                   connectDots=F,
                                   xlab='time steps', ylab='y',
                                   pch=1,
                                   col=c('black','grey')){
  # x: predictor matrix (nXp), 
  #   n is number of observations
  #   p is number of predictors
  #   ng is number of gibbs sampling
  # y: observation time series (nX1)
  # beta: beta coefficients array of (ngxp)
  # sigma process error array of (ngx1)
  # nTrends: number of sampling plots when plotType='all
  
  # y(t) = y(t-1) + x(t-1)*beta + process.error 
  #  process.error ~ N(0, sigma)
  #  observation.error ~ N(0, tau)
  wCon <- rep(0, length(y))
  wCon[startPoints] <- 1
  wCon <- cumsum(wCon)
  
  n <- length(y)
  p <- ncol(x)
  pieces <- cbind(startPoints, c(startPoints[-1]-1, n))
  
  yPred <- matrix(NA, length(y), nTrends)
  for(trend in 1:nTrends){
    beta.Sample <- as.numeric(beta[sample(1:nrow(beta),1 ),])
    sigma.Sample <- as.numeric(sigma[sample(1:length(sigma),1 )])
    
    dy <- rnorm(mean = x%*%beta.Sample, sd = sigma.Sample, n)
    dy <- c(y[1], dy[-n])
    
    yPred[,trend] <- pieceWiseCumSum(dy, startPoints)
    if(!is.null(startPoints))yPred[startPoints,trend] <- y[startPoints]
  }
  if(plotFlag)
  {
    #t <- 1:n
    if(is.null(yLim)) yLim <- range(yPred, na.rm = T)
    plot(t, y, ylim= yLim, xlim=xLim, xlab=xlab, ylab=ylab, col=col[1], pch=pch)
    for(trend in 1:nTrends){
      #lines(t, yPred[,trend], col=col[2])
      for(i in 1:nrow(pieces)){
        st <- pieces[i,1]
        en <- pieces[i,2]
        lines(t[st:en], yPred[st:en, trend], col=col[2])
      }
    }
    yPredQuant <- t(apply(yPred, MARGIN = 1, FUN = quantile, probs=c(.025,.50,.975)))
    lines(t[st:en], yPredQuant[st:en,1], lwd=1.5, lty=2)
    lines(t[st:en], yPredQuant[st:en,2], lwd=1.5, lty=1)
    lines(t[st:en], yPredQuant[st:en,3], lwd=1.5, lty=2)
    if(connectDots)lines(t, y, type = 'l', col=col[2])
    points(t, y, col=col[1], pch=pch)
    #if(!is.null(startPoints))points(startPoints, y[startPoints], col=col[3], pch =19)
  }
  
  invisible(yPred)
}


stateSpacePlot <- function(y, yGibbs,startPoints = 1, title.txt='',xlab='time steps', ylab='y',
                           wMissing=NULL, yLim=NULL, col=c('black','red','grey','blue')){
  n <- length(y)
  pieces <- cbind(startPoints, c(startPoints[-1]-1, length(y)))
  confInterval <- apply(yGibbs, 2, quantile, probs=c(.025,.5,.975))
  
  cols <- rep(col[1], length(y))
  if(!is.null(wMissing))cols[wMissing] <- col[2]
  if(is.null(yLim)) yLim <- range(confInterval)
  
  plot(y, col=cols, ylim=yLim, main = title.txt, xlab=xlab, ylab=ylab)
  #lines(confInterval[2,])
  #polygon(c(1:n, n:1), c(confInterval[1,], rev(confInterval[3,])), col=paste0(col2hex(col[3]),'80'), border = F)
  for(i in 1:nrow(pieces)){
    st <- pieces[i,1]
    en <- pieces[i,2]
    lines(st:en, confInterval[2,st:en])
    polygon(c(st:en, en:st), c(confInterval[1,st:en], confInterval[3,en:st]), col=paste0(col2hex(col[3]),'80'), border = F)
  }
  #lines(confInterval[1,], col=col[3])
  #lines(confInterval[3,], col=col[3])
  segments(x0 = 1:n, y0 = confInterval[1,],
           x1 = 1:n, y1 = confInterval[3,],
           col = col[4])
  points(y, col=cols)
}