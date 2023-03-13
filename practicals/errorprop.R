library(geostats)

LL <- function(lambda,k){
  k * log(lambda) - lambda - sum(1:k)
}

N <- 100
lam <- seq(from=0,to=20,length.out=N)
loglik <- rep(0,N)
for (i in 1:N){
  loglik[i] <- LL(lambda=lam[i],k=4)
}
plot(lam,loglik,type='l',xlab=expression(lambda),ylab='LL')

o <- optim(par=1,f=LL,k=4,control=list(fnscale=-1),method='BFGS')