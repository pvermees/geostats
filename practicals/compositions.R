library(geostats)

uv <- alr(ACNK)
plot(uv)
mu <- colMeans(uv)
covmat <- cov(uv)
points(x=mu[1],y=mu[2],pch=22,bg='black')
ell <- ellipse(mu,covmat)
polygon(ell)

ternary(ACNK,labels=c(expression('Al'[2]*'O'[3]),
                      expression('CaO+Na'[2]*'O'),
                      expression('K'[2]*'O')))

ternary(alr(mu,inverse=TRUE),add=TRUE,type='p',pch=22,bg='black')
ternary(alr(ell,inverse=TRUE),add=TRUE,type='l')

comp <- clr(major)
pc <- prcomp(comp)
biplot(pc)

library(MASS)
affinity <- FAM[,1]
comp <- alr(FAM[,-1])
ld <- lda(x=comp,grouping=affinity)
newrock <- data.frame(F=1,A=8,M=0.1)
newcomp <- alr(newrock)
pr <- predict(object=ld,newdata=newcomp)