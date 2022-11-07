library(geostats)

op <- par(mfrow=c(1,4)) # op = old parameter values
plot(anscombe$x1,anscombe$y1)
plot(anscombe$x2,anscombe$y2)
plot(anscombe$x3,anscombe$y3)
plot(anscombe$x4,anscombe$y4)
par(op)                 # restore old parameter values

np <- 4 # np = 'number of panels'
op <- par(mfrow=c(1,np)) # op = 'old parameters'
for (i in 1:np){
    plot(anscombe[,i],anscombe[,i+np])
}

titles <- c('I','II','III','IV')
for (i in 1:np){
plot(anscombe[,i],anscombe[,i+np],
xlab='x',ylab='y',pch=19,main=titles[i])
}
par(op)

data(catchments,package='geostats')

counts <- table(catchments$lithology)
barplot(counts)

pH <- catchments$pH
hist(pH)
rug(pH)

op <- par(mfrow=c(1,2))
hist(pH,breaks=5)
hist(pH,breaks=10)
par(op)

hist(pH,breaks=seq(from=3,to=7,by=0.5))
hist(pH,breaks=seq(from=3.25,to=6.75,by=0.5))

dens <- density(pH)
plot(dens)
rug(pH)

lc <- log(catchments$CaMg)
d <- density(lc)
plot(d)

veg <- catchments$vegetation/100
lv <- logit(veg)
ld <- density(lv)
plot(ld)

d <- logit(ld,inverse=TRUE)
plot(d)

x <- faithful[,'waiting']
y <- faithful[,'eruptions']
op <- par(mfrow=c(2,1))
plot(density(x),xlab='minutes',main='waiting time')
rug(x)
plot(density(y),xlab='minutes',main='eruption duration')
rug(y)
par(op)

library(MASS)
kde2 <- kde2d(x,y)
contour(kde2)
points(x,y)

cdf <- ecdf(pH)
plot(cdf)

plot(cdf,verticals=TRUE,pch=NA)

cdf(4.5)
