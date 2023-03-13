library(geostats)

d <- rnorm(n=100,mean=50,sd=5)
hist(d)

library(MASS)
m <- c(10,20)
s <- matrix(data=c(2,-3,-3,6),nrow=2,ncol=2)
xy <- mvrnorm(n=200,mu=m,Sigma=s)
plot(xy)

op <- par(mfrow=c(1,2))
m <- 50
s <- 5
x <- seq(from=25,to=75,length.out=100)
f <- dnorm(x=x,mean=m,sd=s)
plot(x=x,y=f,type='l',main='PDF')
P <- pnorm(q=x,mean=m,sd=s)
plot(x=x,y=P,type='l',ylab='P(X<x)',main='CDF')
par(op)