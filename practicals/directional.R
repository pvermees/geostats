library(geostats)

circle.plot(striations,degrees=TRUE)
circle.points(mean(striations),pch=22,
              bg='white',degrees=TRUE)
md <- meanangle(striations,degrees=TRUE)
circle.points(md,pch=21,bg='black',degrees=TRUE)

R <- Rbar(striations,degrees=TRUE)
k <- Rbar2kappa(R)
a <- seq(from=0,to=360,length.out=100)
d <- vonMises(a,mu=md,kappa=k,degrees=TRUE)
plot(a,d,type='l')
rug(striations)

stereonet(trd=palaeomag$decl,plg=palaeomag$incl,
          option=1,degrees=TRUE,show.grid=FALSE)

stereonet(trd=fault$strike,plg=fault$dip,
          option=2,degrees=TRUE,show.grid=FALSE)

stereonet(trd=30,plg=-20,degrees=TRUE,
          option=3,show.grid=TRUE,wulff=FALSE)

stereonet(trd=c(0,0,90,180,270),plg=c(90,10,10,10,10),
          coneAngle=rep(10,5),option=4,degrees=TRUE,
          add=TRUE,wulff=FALSE)

A <- palaeomag$decl
D <- palaeomag$incl
meanpalaeomag <- meanangle(trd=A,plg=D,option=1,degrees=TRUE)
stereonet(trd=A,plg=D,degrees=TRUE,
          show.grid=FALSE,pch=16,col='grey')
stereonet(trd=meanpalaeomag[1],plg=meanpalaeomag[2],
          option=1,degrees=TRUE,add=TRUE,pch=16,col='black')

S <- fault$strike
D <- fault$dip
meanfault <- meanangle(trd=S,plg=D,option=2,degrees=TRUE)
stereonet(trd=S,plg=D,option=2,degrees=TRUE,
          show.grid=FALSE,pch=16,col='grey')
stereonet(trd=meanfault[1],plg=meanfault[2],
          option=2,degrees=TRUE,add=TRUE,pch=16,col='black')

circle.plot(pebbles,degrees=TRUE)
m <- meanangle(pebbles,degrees=TRUE,orientation=TRUE)
circle.points(m,degrees=TRUE,pch=16,cex=2)