# stateSpaceTemporalPlot <- function(x, y, beta, 
#                                    sigma=0, 
#                                    plotType = 'all',
#                                    nPlot=10){
  # x: predictor matrix (nXp), 
  #   n is number of observations
  #   p is number of predictors
  #   ng is number of gibbs sampling
  # y: observation time series (nX1)
  # beta: beta coefficients array of (ngxp)
  # sigma process error array of (ngx1)
  # plotType: factor of "all": start from T=0 and predict the whole wide range
  #                     "1to1": one time step ahead, T=t based on T=t-1
  # nPlot: number of sampling plots when plotType='all
  
  # y(t) = y(t-1) + x(t-1)*beta + process.error 
  #  process.error ~ N(0, sigma)
  #  observation.error ~ N(0, tau)
  
  n <- length(z)
  p <- ncol(x)
  
  
  for(trend in 1:nPlot){
    y0 <- y[1]
    beta.Sample <- beta[sample(1:nrow(beta),1 ),]
    sigma.Sample <- sigma[sample(1:length(sigma),1 )]
    
    
  }
# }