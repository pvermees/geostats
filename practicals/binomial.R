library(geostats)

d <- rbinom(n=50,size=10,prob=0.5)
hist(d)

k <- 0:10
pmf <- dbinom(x=k,size=10,prob=0.5)
barplot(height=pmf,names.arg=k)

cdf <- pbinom(q=0:10,size=10,prob=0.5)
plot(cdf,type='s')