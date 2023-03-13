library(geostats)

d <- rpois(n=100,lambda=3.5)
hist(d)

k <- 0:15
pmf <- dpois(x=k,lambda=3.5)
barplot(height=pmf,names.arg=k)

cdf <- ppois(q=k,lambda=3.5)
plot(cdf,type='s')