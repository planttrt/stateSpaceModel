library('truncnorm')
N <- 100
b <- c(.1, .4)
sig <- .1
sig2 <- .01

p <- length(b)
x <- cbind(1, matrix(runif((p-1)*N),N, p-1))
dy <- rbind(0, x[-N,]%*%b) + rtruncnorm(N, 0, Inf, 0, sig)
dy <- rbind(0, x[-N,]%*%b) + rnorm(N, 0, sig)
y <- cumsum(dy)
z <- rnorm(length(y), y, sig2)
plot(y)



library('rjags')

jags <- jags.model('modelSS.bugs',
                   data = list('x' = x,
                               'z' = z,
                               'N' = N,
                               'p'= ncol(x),
                               connectFore=1,
                               connectBoth=2:(N-1),
                               connectBack=N),
                   n.chains = 4,
                   n.adapt = 100)

update(jags, 1000)

tmp <- jags.samples(jags,
                    c('beta', 'sigma', 'sigma2'),
                    500)

print(tmp)
par(mfrow=c(3,1))
plot(tmp$beta)  
abline(h=b, col='red')

plot(tmp$sigma, type = 'l')
abline(h=sig, col='red')

plot(tmp$sigma2)
abline(h=sig2, col='red')
