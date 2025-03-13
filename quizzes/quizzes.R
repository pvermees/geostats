graphics.off()

setwd('~/git/geostats/quizzes/')
source('../build/notes/helper.R')

install.packages('~/git/geostats/package',repos=NULL,type='source')

library(geostats)

pars <- function(mar=c(2.5,2.3,0.5,0.25),mgp=c(1.5,0.5,0),mfrow=c(1,1)){
    par(list(mar=mar,mgp=mgp,mfrow=mfrow))
}

cairo <- function(file,width,height,family="serif",pointsize=13,...){
    cairo_pdf(file=file,width=width,height=height,
              family=family,pointsize=pointsize,...)
}

options(warn=0)

# PDFs and CDFs for quizzes:
plotpdf <- function(x,d,nr=500,bw="nrd0",xlim=c(0,1)){
    Xd <- rdata(x,d,nr)
    dens <- density(Xd$X,bw=bw)
    xlim[1] <- min(xlim[1],dens$x[1])
    xlim[2] <- max(xlim[2],tail(dens$x,1))
    plot(x=c(dens$x[1],x[1],x,tail(x,1),tail(dens$x,1)),
         y=c(0,0,Xd$dscaled,0,0),type='l',col='grey50',lwd=3,
         bty='n',xlab='x',ylab='Density',
         xlim=xlim,ylim=c(0,max(dens$y,Xd$dscaled)))
    lines(dens$x,dens$y)
    rug(Xd$X)
}
plotcdf <- function(x,d,nr=500){
    Xd <- rdata(x,d,nr)
    dx <- diff(x)
    foo <- cumsum(d*c(dx[1],dx))
    cdf <- foo/tail(foo,n=1)
    plot(x,cdf,type='l',bty='n',xlab='x',ylab='P(x<X)',col='grey50',lwd=3)
    plot(ecdf(Xd$X),verticals=TRUE,pch=NA,add=TRUE,main='')
    rug(Xd$X)
}
rdata <- function(x,d,nr){
    if (length(unique(d))==2){ # delta function
        ud <- unique(d)[2]
        i <- which(d==ud)[1]
        X <- rep(x[i],nr)
        dscaled <- 20*d/max(d)
    } else {
        dx <- c(x[2]-x[1],diff(x))
        area <- sum(d*dx)
        dscaled <- d/area
        cdf <- cumsum(dx*dscaled)
        set.seed(1)
        y <- runif(nr)
        dat <- approx(x=c(0,cdf,1),y=c(x[1],x,tail(x,1)),xout=y)
        X <- dat$y
    }
    list(X=X,dscaled=dscaled)
}
quiz1question <- function(qn=1){
    dx <- 1e-4
    pd <- x <- bw <- list()
    x[[1]] <- seq(from=0,to=1,length.out=100)
    pd[[1]] <- rep(1/100,100)
    bw[[1]] <- "nrd0"
    x[[2]] <- c(seq(from=0,to=0.5-dx,length.out=50),
                0.5,
                seq(from=0.5+dx,to=1,length.out=50))
    pd[[2]] <- c(rep(0,50),1,rep(0,50))
    bw[[2]] <- 0.02
    x[[3]] <- seq(from=0,to=1,length.out=100)
    pd[[3]] <- 0.3*dnorm(x[[3]],mean=0.3,sd=0.05) +
        0.7* dnorm(x[[3]],mean=0.7,sd=0.05)
    bw[[3]] <- "nrd0"
    x[[4]] <- x[[3]]
    pd[[4]] <- 0.7*dnorm(x[[3]],mean=0.3,sd=0.05) +
        0.3* dnorm(x[[3]],mean=0.7,sd=0.05)
    bw[[4]] <- "nrd0"
    x[[5]] <- seq(from=0,to=1,length.out=100)
    pd[[5]] <- x[[5]]/sum(x[[5]])
    bw[[5]] <- "nrd0"
    if (TRUE){ # Moodle
        pars(mfrow=c(5,2),mar=c(3,3,0.5,1))
        plotpdf(x[[qn]],pd[[qn]],bw=bw[[qn]],xlim=c(0,1))
        legend('topleft',legend='1)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plotcdf(x[[3]],pd[[3]])
        legend('topleft',legend='a)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plot.new()
        plotcdf(x[[1]],pd[[1]])
        legend('topleft',legend='b)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plot.new()
        plotcdf(x[[5]],pd[[5]])
        legend('topleft',legend='c)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plot.new()
        plotcdf(x[[2]],pd[[2]])
        legend('topleft',legend='d)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plot.new()
        plotcdf(x[[4]],pd[[4]])
        legend('topleft',legend='e)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
    } else { # Slides
        pars(mfrow=c(3,3),mar=c(3,3,0.5,1))
        plotpdf(x[[qn]],pd[[qn]],bw=bw[[qn]],xlim=c(0,1))
        legend('topleft',legend='1)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plotcdf(x[[3]],pd[[3]])
        legend('topleft',legend='a)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plotcdf(x[[1]],pd[[1]])
        legend('topleft',legend='b)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plot.new()
        plotcdf(x[[5]],pd[[5]])
        legend('topleft',legend='c)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plotcdf(x[[2]],pd[[2]])
        legend('topleft',legend='d)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plot.new()
        plotcdf(x[[4]],pd[[4]])
        legend('topleft',legend='e)',bty='n',adj=c(1,0),xpd=NA,cex=1.2)
        plot.new()
    }
}

for (i in 1:5){
    fn <- paste0('q1q',i,'.png')
    if (TRUE){ # Moodle
        png(filename=fn,width=500,height=800,pointsize=14)
    } else { # Slides
        png(filename=fn,width=1400,height=1000,pointsize=24)
    }
    quiz1question(qn=i)
    dev.off()
}

png(filename='boxplot.png',width=400,height=200,pointsize=14)
pars(mar=c(2,0.5,0.5,0.5))
dat <- c(seq(from=0,to=50,length.out=25),
         seq(from=51,to=100,length.out=25),
         seq(from=101,to=200,length.out=25),
         seq(from=201,to=400,length.out=25)
         )
boxplot(dat,horizontal=TRUE,xaxt='n')
axis(1,at=c(0,50,100,200,400))
dev.off()

png(filename='ecdf.png',width=400,height=400,pointsize=14)
pars(mar=c(3,3,0.5,0.5))
plot(ecdf(1:9),verticals=TRUE,pch=NA,main=NA,xaxt='n')
axis(1,at=1:9)
dev.off()

if (FALSE){
    tab <- rbind(dbinom(x=0:10,size=10,prob=0.4),
                 pbinom(q=0:10,size=10,prob=0.4))
    colnames(tab) <- 0:10
    signif(tab,2)
}
if (FALSE){
    kk <- 60
    nn <- 100
    qq <- c(0.01,0.025,0.05,0.1,0.5,0.9,0.95,0.975,0.99)
    pp <- c(binom.test(kk,nn,conf.level=0.99,alternative='two.sided')$conf.int[1],
            binom.test(kk,nn,conf.level=0.95,alternative='two.sided')$conf.int[1],
            binom.test(kk,nn,conf.level=0.95,alternative='greater')$conf.int[1],
            binom.test(kk,nn,conf.level=0.95,alternative='less')$conf.int[2],
            binom.test(kk,nn,conf.level=0.95,alternative='two.sided')$conf.int[2],
            binom.test(kk,nn,conf.level=0.99,alternative='two.sided')$conf.int[2]
            )
    out <- NULL
    for (p in pp){
        out <- rbind(out,qbinom(qq,prob=p,size=nn))
    }
    colnames(out) <- qq
    rownames(out) <- signif(pp,3)
}

png(filename='dinosaurs2sided.png',width=1000,height=500)
pars(mfrow=c(1,2))
binomhist(nn=50,kk=18,H0=0.5,nsides=2,xlab='x = # male dinosaurs')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=50,kk=18,H0=0.5,nsides=2,xlab='x = # male dinosaurs')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

png(filename='dinosaursless.png',width=1000,height=500)
pars(mfrow=c(1,2))
binomhist(nn=50,kk=18,H0=0.5,nsides=1,xlab='x = # male dinosaurs')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=50,kk=18,H0=0.5,nsides=1,xlab='x = # male dinosaurs')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

png(filename='dinosaursmore.png',width=1000,height=500)
pars(mfrow=c(1,2))
binomhist(nn=50,kk=32,H0=0.5,nsides=-1,xlab='x = # female dinosaurs')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=50,kk=32,H0=0.5,nsides=-1,xlab='x = # female dinosaurs')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

png(filename='ci.png',width=400,height=400,pointsize=14)
pars(mar=c(3,3,0,0))
xx <- c(1,2,3,4)
yy <- c(20,22,19,21)
dy <- c(4,1,3,2)
plot(xx,yy,type='n',ann=FALSE,bty='n',xaxt='n',
     xlim=c(0.5,4.5),ylim=range(c(yy-dy,yy+dy)))
arrows(xx,yy-dy,xx,yy+dy,code=3,angle=90)
axis(1,at=1:4,labels=c('a','b','c','d'))
dev.off()

png(filename='QQ.png',width=800,height=800,pointsize=14)
data(DZ,package='geostats')
ns <- 5
p <- par(mfrow=c(ns,ns),mar=c(3,3,0.5,0.5),mgp=c(1.5,0.5,0))
for (i in 1:ns){
    for (j in 1:ns){
        if (i==j){
            plot.new()
        } else {
            qqplot(DZ[[i]],DZ[[j]],xlab=paste('sample',i),ylab=paste('sample',j))
        }
    }
}
par(p)
dev.off()

png(filename='pt.png',width=400,height=300,pointsize=14)
p <- par(mar=c(3,3,0.5,0.5),mgp=c(1.5,0.5,0))
x <- seq(from=-5,to=5,length.out=100)
y <- pt(x,df=3)
plot(x,y,type='l',xlab='t',ylab='F')
lines(range(x),rep(0.025,2),lty=2)
lines(range(x),rep(0.975,2),lty=2)
lines(rep(qt(0.025,df=3),2),c(0,1),lty=3)
lines(rep(qt(0.975,df=3),2),c(0,1),lty=3)
par(p)
dev.off()

png(filename='majorqq.png',width=1200,height=400,pointsize=14)
pars(mfrow=c(2,5),mar=c(3,3,1,0.5))
data(major,package='geostats')
for (cname in colnames(major)){
    qqnorm(major[,cname],main=cname)
    qqline(major[,cname])
}
dev.off()

if (FALSE){
    tab <- matrix(0,7,10)
    for (i in 1:10){
        tab[,i] <- signif(qt(c(0.01,0.025,0.05,0.5,0.95,0.975,0.99),df=i),3)
    }
}

png(filename='comparingdistributions-2-3.png',
    width=400,height=400,pointsize=14)
pars()
chi2 <- seq(from=0,to=10,length.out=100)
plot(x=chi2,y=pchisq(chi2,df=4),type='l',
     xlab=expression(chi^2),ylab=('F'))
lines(x=range(chi2),y=rep(0.95,2),lty=3)
lines(x=rep(qchisq(0.95,df=4),2),y=c(0,1),lty=3)
dev.off()

png(filename='regression-1-2.png',
    width=400,height=400,pointsize=14)
E <- rbind(c(1,0.1),c(0.1,1))
e <- MASS::mvrnorm(500,rep(0,2),E)
colnames(e) <- c('x','y')
fit <- lm(e[,'y']~e[,'x'])
tit <- paste0('r=',signif(sqrt(summary(fit)$r.squared),2),
              ', p-value=',signif(summary(fit)$coefficients[2,4],2))
plot(e,pch=16,bty='n',main=tit)
dev.off()

png(filename='regression-1-3a.png',
    width=400,height=400,pointsize=14)
par(mfrow=c(2,2),mar=c(3,3,0.5,0.5),mgp=c(1.5,0.5,0))
E <- rbind(c(1,-0.7),c(-0.7,1))
e <- MASS::mvrnorm(50,rep(0,2),E)
x <- 10 + e[,1]
y <- 10 + e[,2]
z <- rnorm(50,mean=10,sd=3)
plot(x,y,pch=16)
plot(x,z,pch=16)
plot(y,z,pch=16)
plot(x/z,y/z,pch=16)
dev.off()

png(filename='regression-1-3b.png',
    width=400,height=400,pointsize=14)
par(mfrow=c(2,2),mar=c(3,3,0.5,0.5),mgp=c(1.5,0.5,0))
E <- rbind(c(1,-0.95),c(-0.95,1))
e <- MASS::mvrnorm(50,rep(0,2),E)
x <- 10 + e[,1]
y <- 10 + e[,2]
z <- rnorm(50,mean=50,sd=1)
plot(x,y,pch=16)
plot(x,z,pch=16)
plot(y,z,pch=16)
plot(x/z,y/z,pch=16)
dev.off()

png(filename='unsupervised-pca.png',
    width=400,height=400,pointsize=14)
pars(mar=c(3,3,2,2))
dat <- read.csv('Zanzibar.csv',header=TRUE)
pc <- prcomp(dat[,-c(1,2)],scale=TRUE)
biplot(pc,xlabs=dat$channel,xlim=c(-0.5,0.7))
dev.off()

set.seed(5)
#png(filename='unsupervised-hclust.png',
#    width=250,height=750,pointsize=14)
#pars(mfrow=c(4,1))
pdf(file='../slides/unsupervised-hclust.pdf',
    width=6,height=6)
pars(mfrow=c(2,2))
labs <- c('(a)','(b)','(c)','(d)')
xy <- list()
d <- list()
for (i in 1:4){
    xy[[i]] <- cbind(runif(4,min=0,max=100),
                     runif(4,min=50,max=200))
    d[[i]] <- round(dist(xy[[i]]))
    hc <- hclust(d[[i]])
    plot(hc,main='',sub='',xlab='')
    legend('topright',labs[i],xpd=NA,bty='n',cex=1.2)
}
dev.off()
print(d[[1]])
pdf(file='unsupervised-hclust.pdf')
plot(xy[[1]],type='n',xlab='x',ylab='y')
text(xy[[1]],labels=1:4)
dev.off()

pdf(file='crimeclust.pdf')
crime <- scale(USArrests)
fit <- hclust(dist(crime))
plot(fit,xlab='',sub='')
dev.off()

if (FALSE){
    fit <- kmeans(crime,centers=4)
    fit$cluster[c('Alaska','Texas')]

    data(training,package='geostats')
    my.control <- rpart.control(xval=10, cp=0, minsplit=1)
    tree <- rpart::rpart(affinity ~ ., data=training, control=my.control)
    rpart::printcp(tree)
    tree <- rpart::rpart(affinity ~ ., data=training)
    rpart::printcp(tree)
}

png(filename='compositional-biplot.png',
    width=400,height=400,pointsize=14)
pars(mar=c(3,3,2,2))
set.seed(4)
lrcomp <- MASS::mvrnorm(5,mu=c(0,0,0),Sigma=rbind(c(2,1,2),c(1,3,1),c(2,1,5)))
comp <- round(100*clr(lrcomp,inverse=TRUE))
lrcomp2 <- clr(comp)
colnames(lrcomp2) <- c('Ca','Mg','Fe')
pc <- prcomp(lrcomp2)
biplot(pc)
dev.off()

png(filename='wulff.png',
    width=400,height=400,pointsize=14)
pars(mar=rep(1.3,4))
stereonet(wulff=TRUE)
dev.off()

png(filename='trdplgsolution.png',
    width=600,height=300,pointsize=14)
pars(mar=rep(1.3,4),mfrow=c(1,2))
stereonet(trd=30,plg=10,degrees=TRUE,wulff=FALSE,
          show.grid=FALSE,pch=21,bg='black')
legend('topleft',legend='a)',bty='n',adj=c(2,0))
stereonet(trd=30,plg=10,degrees=TRUE,wulff=TRUE,
          show.grid=FALSE,pch=21,bg='black')
legend('topleft',legend='b)',bty='n',adj=c(2,0))
dev.off()

pdf(file='trdplgsolution.pdf',width=6,height=3)
pars(mar=rep(1.3,4),mfrow=c(1,2))
stereonet(trd=30,plg=10,degrees=TRUE,wulff=FALSE,
          show.grid=TRUE,pch=21,bg='black')
legend('topleft',legend='a)',bty='n',adj=c(2,0))
stereonet(trd=30,plg=10,degrees=TRUE,wulff=TRUE,
          show.grid=TRUE,pch=21,bg='black')
legend('topleft',legend='b)',bty='n',adj=c(2,0))
dev.off()

png(filename='meuse.png',width=600,height=600,pointsize=14)
data('meuse',package='geostats')
colourplot(X=meuse$x,Y=meuse$y,Z=log(meuse$zinc),
           plot.title=title(main='Meuse',xlab='Easting',ylab='Northing'),
           key.title=title(main='log(Zn)'))
points(178800,330750,pch=21,bg='white',cex=4)
text(178800,330750,'a',cex=2)
points(179450,331700,pch=21,bg='white',cex=4)
text(179450,331700,'b',cex=2)
points(179500,332500,pch=21,bg='white',cex=4)
text(179500,332500,'c',cex=2)
dev.off()
