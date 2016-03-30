

myrmvnorm <- function (nn, mu, sigma){
  
  
  testv <- try(svd(sigma),T)
  
  if( inherits(testv,'try-error') ){
    
    ev <- eigen(sigma, symmetric = TRUE)
    testv <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
  } else {
    testv <- t(testv$v %*% (t(testv$u) * sqrt(testv$d)))
  }
  
  retval <- matrix(rnorm(nn * ncol(sigma)), nn) %*% testv
  retval + mu
}


invMat <- function(SS,NEARPD=F){  #matrix inversion, if NEARPD find closest PD matrix
  
  require(Matrix)
  
  testv <- try(chol(SS),T)
  
  if( inherits(testv,'try-error') ){
    message('near pos definite used in invMat')
    if(NEARPD){
      require(Matrix)
      SS     <- as.matrix( nearPD(SS)$mat )
      testv <- try(chol(SS),T)
    }
  }
  
  chol2inv(testv)
}
bUpdateNorm <- function(xx,yy,b,
                        priorB=matrix(ncol(xx)*0,ncol(xx)),priorIVB=diag(1/100,ncol(xx)),
                        loB=NULL,hiB=NULL,sigma){
  
  V <- invMat( crossprod(xx)/sigma + priorIVB )
  v <- crossprod(xx,yy)/sigma + priorIVB%*%priorB
  if( is.null(loB) & is.null(hiB) )return( t( myrmvnorm(1,t(V%*%v),V) ) )
  
  tnorm.mvt(V%*%v,V%*%v,V,loB,hiB)
}


updateVariance <- function(yy,mu,s1=1,s2=1,lo=NULL,hi=NULL){
  
  # if yy is a matrix, one variance for each column
  
  require(pscl)
  
  tiny <- 1e-10
  
  k <- 1
  res <- (yy - mu)^2
  nn  <- length(yy)
  
  if(is.matrix(yy)){
    k  <- ncol(yy)
    nn <- nrow(yy)
    sr <- colSums(res)
  }else{
    sr <- sum(res)
  }
  
  
  u1 <- s1 + nn/2
  u2 <- s2 + .5*sr
  
  if(is.null(lo) & is.null(hi))return( 1/rgamma(k,u1,u2) )
  
  if(is.null(lo))lo <- 0
  
  rtrunc(k,lo,hi,u1,u2,'igamma')
  
}

tnorm <- function(n,lo,hi,mu,sig){   
  
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == Inf]  <- lo[z == Inf]
  z[z == -Inf] <- hi[z == -Inf]
  z
}


predVsObs <- function(o,p,ylim=range(p,na.rm=T),xlab=' ',ylab=' ',colors=rep(1,length(o))){ 
  
  #o  - length n vector of obs or true values
  #p - ng by n matrix of estimates
  
  n <- length(o)
  y <- apply(p,2,quantile,c(.5,.025,.975))
  
  plot(o,y[1,],ylim=ylim,xlab=xlab,ylab=ylab,col=colors)
  for(j in 1:n)lines(c(o[j],o[j]),y[2:3,j],col=colors[j])
  abline(0,1,lty=2)
  invisible(y)
}
