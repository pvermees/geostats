library(geostats)

plot(x=rbsr$RbSr,y=rbsr$SrSr)

r <- cor(x=rbsr$RbSr,y=rbsr$SrSr)

n <- nrow(rbsr)
tstat <- r*sqrt(n-2)/sqrt(1-r^2)

fit <- lm(SrSr ~ RbSr, data=rbsr)
abline(fit)
x <- seq(from=min(rbsr$RbSr),to=max(rbsr$RbSr),length.out=20)
pred <- predict(fit,newdata=data.frame(RbSr=x),
interval="confidence",level=0.95)
matlines(x,pred,lty=1,col='black')
pred <- predict(fit,newdata=data.frame(RbSr=x),
                interval="prediction",level=0.95)
matlines(x,pred,lty=2,col='black')

yfit <- york(x=rbsr)