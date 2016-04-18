library('truncnorm')
source('~/Projects/procVisData/bayesianFunctions.R')
library(rjags)

setConnectBackFore <- function(xyData){
  xyData.tmp <- xyData
  xyData.tmp[,ConnectBack:=0:(.N-1)]
  xyData.tmp[,ConnectFore:=2:(.N+1)]
  
  xyData.tmp$ConnectBack[xyData$Head] <- NA
  xyData.tmp$ConnectFore[xyData$Tail] <- NA
  xyData.tmp$ConnectBack
  
  xyData.tmp
}

pLatentStates.SS <- function(x, z, yg, bg, sg, tg, connect, wNA, TRUNC){ 
  p.obs  <- dnorm(z, yg, tg, log=T)
  
  if(TRUNC){
    p.back <- log(dtruncnorm(yg - yg[connect[,1]], 0, Inf, x[connect[,1],]%*%bg, sg))
    p.fore <- log(dtruncnorm(yg[connect[,2]] - yg, 0, Inf, x%*%bg, sg))
  }else{
    p.back <- dnorm(yg- yg[connect[,1]], x[connect[,1],]%*%bg, sg, log=T)
    p.fore <- dnorm(yg[connect[,2]]- yg, x%*%bg, sg, log=T)
  }
  
  p.obs[is.na(p.obs)] <- 0
  p.back[is.na(p.back)] <- 0
  p.fore[is.na(p.fore)] <- 0
  
  p <- p.obs+p.back+p.fore
}


stateSpace <- function(x, z, connect=NULL, HeadTail=NULL,
                       ng = 1000, burnin = floor(0.5 * ng), 
                       MPstep = .01, nAccept =50, 
                       priorB = rep(0, ncol(x)), priorIVB = diag(1/1000, ncol(x)),
                       tauPriorStrength=NULL, tauPriorMean=NULL, storeLatent=F,
                       sigPriorStrength=NULL, sigPriorMean=NULL, TRUNC =F){
  if(!is.null(connect)){
    Head <- which(is.na(connect[,1]))
    Tail <- which(is.na(connect[,2]))
    Body <- which(!rowSums(is.na(connect)))
  }else if(!is.null(HeadTail)){
    Head <- which(HeadTail[,1])
    Tail <- which(HeadTail[,2])
    Body <- which(!(HeadTail[,1]|HeadTail[,2]))
  }else
  {
    stop('One of connect or HeadTail arguments should be defined.')
  }
  
  if(any(is.na(x))) stop('Matrix x is not allowed to contain NAs!')
  
  wNA <- which(is.na(z))
  
  N <- nrow(x)
  p <- ncol(x)
  
  MPstep <- rep(MPstep, N)
  
  bg <- rep(1, p)
  sg <- tg <- 10
  
  if(is.null(tauPriorMean)|is.null(tauPriorStrength)){
    print('no prior for tau!')
    t1 <- .1
    t2 <- .1
  }else{
    print(paste0('a prior on tau = ', tauPriorMean))
    t1 <- N*tauPriorStrength
    t2 <- (tauPriorMean)^2*(t1-1)
  }
  
  if(is.null(sigPriorMean)|is.null(sigPriorStrength)){
    print('no prior for sigma!')
    s1 <- .1
    s2 <- .1
  }else{
    print(paste0('a prior on sigma = ',sigPriorMean))
    s1 <- N*sigPriorStrength
    s2 <- (sigPriorMean)^2*(s1-1)
  }
  
  
  if(storeLatent) ygibbs <- matrix(NA, nrow = ng, ncol = N)
  
  zPredGibbs <- matrix(NA, nrow = ng, ncol = N)
  
  bgibbs <- matrix(0,ng,p)
  colnames(bgibbs) <- colnames(x)
  
  vgibbs <- matrix(0, ng, 2)
  colnames(vgibbs) <- c('sig','tau')
  
  agibbs <- matrix(NA,ng, N)
  accept <- matrix(0,floor(ng/nAccept), N)
  MPStepsgibbs <- matrix(NA,floor(ng/nAccept), N)
  cnt <- 0
  
  yg <- z
  while(any(is.na(yg))){
    yg [is.na(yg)] <- yg[connect[is.na(yg),1]]
    yg [is.na(yg)] <- yg[connect[is.na(yg),2]]
  }
  
  yg.var <- yg.mean <- yg*0
  
  prog <- txtProgressBar(1,ng,1,style = 3, title = 'Gibbs Sampling ...')
  
  for(g in 1:ng){
    yg1 <- yg[!is.na(connect)[,2]]
    yg2 <- yg[!is.na(connect)[,1]]
    x1 <- x[!is.na(connect)[,2],]
    yy <- yg2 - yg1
    xx <- x1
    
    if(TRUNC){
      ff <- x1%*%bg>0
    }else{
      ff <- rep(T, length(yy))
    }    
    bg <- bUpdateNorm(xx = xx[ff,], yy = yy[ff], 
                      b = bg, sigma = sg^2,
                      priorB = priorB, priorIVB = priorIVB)
    
    sg <- sqrt(updateVariance(yg2, yg1 + x1%*%bg, s1, s2))
    tg <- sqrt(updateVariance(z[-wNA], mu=yg[-wNA], t1, t2))
    
    pnow <- pLatentStates.SS(x, z, yg, bg, sg, tg, connect, wNA=wNA, TRUNC = TRUNC)
    
    pyg <- rnorm(N, yg, MPstep)
    
    pnew <- pLatentStates.SS(x, z, pyg, bg, sg, tg, connect, wNA=wNA, TRUNC = TRUNC)
    
    a <- exp(pnew - pnow)
    ran <- runif(N,0,1)
    a.cond <- ran < a
    a.cond[is.na(a.cond)] <- F
    w <- which(a.cond)
    
    yg[w] <- pyg[w]
    
    if(g >= burnin){
      yg.mean.old <- yg.mean
      yg.mean <- yg.mean + yg/(ng-burnin)
    }
    
    bgibbs[g,] <- bg   #save estimates
    vgibbs[g,]  <- c(sg,tg)
    agibbs[g,]  <- a.cond
    
    if(g%%nAccept==0){
      cnt <- cnt + 1
      accept[cnt,] <- colMeans(agibbs[(cnt-1)*nAccept+(1:nAccept),])
      MPstep <- MPstep*(1.1*(accept[cnt,]>0.325) + 0.99*(accept[cnt,]<0.275) + 1*(accept[cnt,]>=0.275&accept[cnt,]<=0.325))
      MPStepsgibbs[cnt,] <- MPstep
    }
    
    if(storeLatent) ygibbs[g,] <- yg
    
    zPredGibbs[g,] <- rnorm(length(yg), yg, tg)
    
    setTxtProgressBar(prog,g)
  }
  
  if(storeLatent) yg.var <- apply(ygibbs[-(1:burnin),], 2, var)
  if(storeLatent) tmp <- list(bgibbs=bgibbs, sgibbs=vgibbs[,1], tgibbs=vgibbs[,2], nBurnin=burnin, mpSteps=MPStepsgibbs, acceptGibbs = accept, zPredGibbs=zPredGibbs, latentMean =yg.mean, latentStd =yg.var, latentGibbs=ygibbs)
  if(!storeLatent) tmp <- list(bgibbs=bgibbs, sgibbs=vgibbs[,1], tgibbs=vgibbs[,2], nBurnin=burnin, mpSteps=MPStepsgibbs, acceptGibbs = accept, zPredGibbs=zPredGibbs, latentMean =yg.mean)
  tmp$chains <- data.frame(beta=tmp$bgibbs,
                          sigma=tmp$sgibbs,
                          tau=tmp$tgibbs)
  
  tmp
}




ssSimulations <- function(nSites=1000, nTSet=c(3:6), p=2, beta =NULL,
                          sig= .1, tau=.01, miss=0,
                          plotFlag = F, TRUNC = F){
  
  nSamples.Site <- sample(nTSet, nSites, replace = T)
  
  if(length(nTSet)==1) nSamples.Site <- sample(c(nTSet,nTSet), nSites, replace = T)
  
  N <- sum(nSamples.Site)
  
  if(!is.null(beta)) {
    p <- length(beta)
  }else {
    beta <- matrix(2*runif(p)-1)
  }
  
  x <- matrix(runif(p*N), ncol=p, nrow=N)
  
  
  y <- rep(0, N)
  
  sampleSiteNo <- c(0, cumsum(nSamples.Site))
  for(i in 1:nSites){
    y[sampleSiteNo[i]+1] <- runif(1, .2, .3)
    
    for(j in 2:nSamples.Site[i])
      #xxx
      if(TRUNC){
        y[sampleSiteNo[i]+j] <- y[sampleSiteNo[i]+j-1] + tnorm(1, 0,Inf, x[sampleSiteNo[i]+j-1,]%*%beta, sig )
      }else{
        y[sampleSiteNo[i]+j] <- y[sampleSiteNo[i]+j-1] + rnorm(1, x[sampleSiteNo[i]+j-1,]%*%beta, sig )
      }
  }
  
  z <- rnorm(N, y, tau)
  
  
  
  connect <- matrix(NA, nrow = N, ncol = 2)
  colnames(connect) <- c('Back','Fore')
  
  for(i in 1:nSites){
    connect[(sampleSiteNo[i]+1):(sampleSiteNo[i+1]-1),2] <- 
      (sampleSiteNo[i]+2):(sampleSiteNo[i+1])
    
    connect[(sampleSiteNo[i]+2):(sampleSiteNo[i+1]),1] <- 
      (sampleSiteNo[i]+1):(sampleSiteNo[i+1]-1)
  }
  
  # beta <- c(-.1, 10)
  # p <- length(beta)
  
  #z1 <- z.[!order[,2]]
  #z2 <- z.[!order[,1]]
  #   x <- x.[!order[,2],]
  #   z <- z.
  wNA <- sample(1:N, floor(miss*N) )
  z[wNA] <- NA
  if (plotFlag)
  {
    plot(z)
    lines(z)
  }
  
  list(x=x, z=z, y= y, connect=connect, 
       tau = tau, sig=sig, beta=beta, 
       startPoints = which(is.na(connect[,1])), 
       n=nSamples.Site, wNA = wNA, TRUNC=TRUNC)
}

stateSpaceJags <- function(x, z, connect=NULL,
                           HeadTail =NULL,
                           nGibbs = 1000, nPost=nGibbs,
                           TRUNC =F,
                           n.chains=4, n.adapt=100,
                           quiet=F, calcLatentGibbs=F
                           
){
  if(!is.null(connect)){
    Head <- which(is.na(connect[,1]))
    Tail <- which(is.na(connect[,2]))
    Body <- which(!rowSums(is.na(connect)))
  }else if(!is.null(HeadTail)){
    Head <- which(HeadTail[,1])
    Tail <- which(HeadTail[,2])
    Body <- which(!(HeadTail[,1]|HeadTail[,2]))
  }else
  {
    stop('One of connect or HeadTail arguments should be defined.')
  }
  
  model <- ifelse(TRUNC, 
                  '~/Projects/stateSpaceModel/modelSS.trunc.bugs',
                  '~/Projects/stateSpaceModel/modelSS.bugs')
  
  # model <- ifelse(TRUNC, modelSS.trunc.bugs,modelSS.bugs)
  
  ssModel <- jags.model(model, 
                        quiet = quiet, 
                        data = list('x' = x,
                                    'z' = z,
                                    'N' = nrow(x),
                                    'p'= ncol(x),
                                    connectHead=Head,
                                    connectBody=Body,
                                    connectTail=Tail),
                        n.chains = n.chains,
                        n.adapt = n.adapt)
  
  update(ssModel, nGibbs)
  
  ssSamples <- jags.samples(ssModel,c('y','beta', 'sigma', 'sigma2'), nPost )
  
  
  
  # print(ssSamples)
  
  
  ssGibbs.jags <- data.frame(beta=t(apply(ssSamples$beta, c(1,2), mean)),
                             sigma=t(apply(ssSamples$sigma, c(1,2), mean)),
                             tau=t(apply(ssSamples$sigma2, c(1,2), mean)))
  latentGibbs <- NULL
  if(calcLatentGibbs) latentGibbs <- t(apply(ssSamples$y, c(1,2), mean))
  
  return(list(model=ssModel, chains=ssGibbs.jags, 
              latentGibbs = latentGibbs, 
              rawsamples = ssSamples))
}


