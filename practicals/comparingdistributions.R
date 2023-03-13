library(geostats)

qqplot(faithful[,'eruptions'],faithful[,'waiting'])

dat <- faithful[,'eruptions']
qqnorm(dat)
qqline(dat)

gold1 <- c(19.09,19.17,19.31,19.07,19.18)
h <- t.test(gold1,mu=19.30,alternative='less')

gold2 <- c(19.30,19.33,19.15,19.32)
h <- t.test(gold1,gold2)

quakesperyear <- countQuakes(declustered,minmag=5.0,from=1917,to=2016)
obs <- table(quakesperyear)

lam <- mean(quakesperyear)
# expected probabilities:
p <- c(ppois(q=1,lambda=lam),                   # lower tail
       dpois(x=2:11,lambda=lam),                # central segment
       ppois(q=11,lambda=lam,lower.tail=FALSE)) # upper tail
names(p) <- c('<2',2:11,'>11')

pred <- p*length(quakesperyear)

O <- c(sum(obs[1:2]),obs[3:8],sum(obs[9:12]))
E <- c(sum(pred[1:2]),pred[3:8],sum(pred[9:12]))

X2 <- sum((O-E)^2/E)
pval <- pchisq(X2,df=length(E)-2,lower.tail=FALSE)

gold1 <- c(19.09,19.17,19.31,19.07,19.18)
gold2 <- c(19.30,19.33,19.15,19.32)
h <- wilcox.test(x=gold1,y=gold2)

river <- DZ[['Y']]
dune <- DZ[['5']]
h <- ks.test(x=river,y=dune)