graphics.off()

setwd('~/Documents/Programming/R/geostats/build/notes/')
source('helper.R')

install.packages('~/Documents/Programming/R/geostats/package',
                 repos=NULL,type='source')

pars <- function(mar=c(2.5,2.3,0.5,0.25),mgp=c(1.5,0.5,0),mfrow=c(1,1)){
    par(list(mar=mar,mgp=mgp,mfrow=mfrow))
}

cairo <- function(file,width,height,family="serif",pointsize=13,...){
    cairo_pdf(file=file,width=width,height=height,
              family=family,pointsize=pointsize,...)
}

options(warn=0)

# anscombe data

cairo(file='../../figures/anscombe.pdf',width=7,height=1.75)
pars(mfrow=c(1,4),mar=c(2.5,2.2,1,0.2))
plot(anscombe$x1,anscombe$y1,xlab='x',ylab='y',pch=19); title('I')
plot(anscombe$x2,anscombe$y2,xlab='x',ylab='y',pch=19); title('II')
plot(anscombe$x3,anscombe$y3,xlab='x',ylab='y',pch=19); title('III')
plot(anscombe$x4,anscombe$y4,xlab='x',ylab='y',pch=19); title('IV')
dev.off()

cairo(file='../../figures/clasts.pdf',width=4,height=3)
pars()
clasts <- c(10,5,6,20)
names(clasts) <- c('granite','basalt','gneiss','quartzite')
barplot(clasts,col='white')
dev.off()

cairo(file='../../figures/binwidth.pdf',width=6,height=3)
pars(mfrow=c(1,2))
pH <- c(6.2,4.4,5.6,5.2,4.5,5.4,4.8,5.9,3.9,3.8,
        5.1,4.1,5.1,5.5,5.1,4.6,5.7,4.6,4.6,5.6)
p <- par(mfrow=c(1,2))
hist(pH,breaks=seq(from=3,to=7,by=1),xlim=c(3,7),main='')
rug(pH)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(1,1))
hist(pH,breaks=seq(from=3,to=7,by=0.5),xlim=c(3,7),main='')
rug(pH)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(1,1))
dev.off()

cairo(file='../../figures/binpos.pdf',width=6,height=3)
pars(mfrow=c(1,2))
hist(pH,breaks=seq(from=3,to=7,by=0.5),xlim=c(3,7),main='')
rug(pH)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(1,1))
hist(pH,breaks=seq(from=3.25,to=6.75,by=0.5),xlim=c(3,7),main='')
rug(pH)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(1,1))
dev.off()

cairo(file='../../figures/rectKDE.pdf',width=3,height=3)
pars()
dat <- c(-0.5,0.3,0.6)
col <- c('gray40','gray60','gray80')
bw <- 1
dens <- density(dat,kernel='rectangular',n=2^13,bw=bw/sqrt(3))
X <- dens$x
Y <- 3*dens$y/max(dens$y)
plot(X,Y,main='',type='n',yaxt='n',ylab='Frequency')
axis(side=2,at=0:3)
for (i in 1:length(dat)){
    x <- dat[i] + c(-1,-1,1,1)*bw
    y <- c(0,1,1,0)
    lines(x,y,lty=2,col=col[i])
}
lines(X,Y)
rug(dat,col=col)
dev.off()

cairo(file='../../figures/pHrectKDE.pdf',width=3,height=3)
pars()
dens <- density(pH,kernel='rectangular',n=2^13)
plot(dens,main='',xlab='pH',zero.line=FALSE)
rug(pH)
dev.off()

cairo(file='../../figures/gaussKDE.pdf',width=3,height=3)
pars()
bw <- 1
dens <- density(dat,n=2^13,bw=bw/sqrt(3))
plot(dens,main='',xlab='X',zero.line=FALSE)
x <- dens$x
for (i in 1:length(dat)){
    y <- dnorm(x,dat[i],bw)/3
    lines(x,y,lty=2,col=col[i])
}
rug(dat,col=col)
dev.off()

cairo(file='../../figures/pHgaussKDE.pdf',width=3,height=3)
pars()
dens <- density(pH)
plot(dens,main='',xlab='pH',zero.line=FALSE)
rug(pH)
dev.off()

cairo(file='../../figures/bandwidth.pdf',width=6,height=3)
pars(mfrow=c(1,2))
dens <- density(pH,bw=0.1)
plot(dens,main='',xlab='pH',zero.line=FALSE)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
rug(pH)
dens <- density(pH,bw=1)
plot(dens,main='',xlab='pH',zero.line=FALSE)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
rug(pH)
dev.off()

cairo(file='../../figures/negativeKDE.pdf',width=3,height=3)
pars()
#clasts <- exp(rnorm(20,1,1))
clasts <- c(0.35, 11.00, 6.00, 1.80, 2.30, 0.59, 8.40, 2.90,
            5.90, 2.10, 1.20, 2.10, 1.10, 1.60, 0.90, 1.70, 3.40, 0.53, 2.20, 7.70)
plot(density(clasts),main='',xlab='clast size (cm)',zero.line=FALSE)
lines(c(0,0),c(0,1),lty=2)
rug(clasts)
dev.off()

cairo(file='../../figures/logKDE.pdf',width=6,height=3)
pars(mfrow=c(1,2))
lclasts <- log(clasts)
ldens <- density(lclasts)
plot(exp(ldens$x),ldens$y,type='l',main='',
     xlab='clast size (cm)',ylab='Density',log='x')
rug(clasts)
legend('topright',legend='a)',bty='n',cex=1.2,adj=c(0,0))
clastdens <- geostats:::exp.density(ldens)
plot(clastdens,main='',xlab='clast size (cm)',
     ylab='Density',xlim=c(0,20),zero.line=FALSE)
rug(clasts)
legend('topright',legend='b)',bty='n',cex=1.2,adj=c(0,0))
dev.off()

cairo(file='../../figures/porosityKDE.pdf',width=3,height=3)
pars()
#u <- c(rnorm(10,-2,1),rnorm(10,2,1))
por <- c(0.058,0.280,0.120,0.270,0.400,0.120,0.038,0.063,0.170,0.160,
         0.950,0.940,0.920,0.880,0.880,0.700,0.920,0.720,0.740,0.840)
u <- log(por/(1-por))
#por <- exp(u)/(exp(u)+1)
plot(density(por),main='',xlab='porosity (fraction)',zero.line=FALSE)
rug(por)
lines(c(0,0),c(0,1),lty=2)
lines(c(1,1),c(0,1),lty=2)
dev.off()

cairo(file='../../figures/logitKDE.pdf',width=6,height=3.1)
pars(mfrow=c(1,2),mar=c(2.5,2.3,2.5,0.2))
ldens <- density(u)
plot(ldens,main='',xlab='porosity (fraction)',xaxt='n',zero.line=FALSE)
rug(u)
ticks <- c(0.01,0.1,0.5,0.9,0.99)
axis(side=1,at=log(ticks/(1-ticks)),labels=ticks)
axis(side=3)
mtext(text='logit[porosity]',side=3,line=1.5)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
pordens <- geostats:::logit.density(ldens,inverse=TRUE)
plot(pordens,main='',xlab='porosity (fraction)',zero.line=FALSE)
rug(por)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(0.5,0))
dev.off()

cairo(file='../../figures/KDE2D.pdf',width=5,height=5)
pars(mar=rep(0,4))
m <- rbind(c(1,2,2,3),c(1,4,5,3),c(1,6,7,3),c(1,8,8,3))
layout(m,widths=c(0.09,0.65,0.25,0.01),heights=c(0.01,0.25,0.65,0.09))
plot.new()
plot.new()
plot.new()
xlim <- c(40,100)
ylim <- c(1.1,5.5)
wdens <- density(faithful[,'waiting'])
plot(wdens,main='',xaxt='n',xlab='',xlim=xlim,zero.line=FALSE)
rug(jitter(faithful[,'waiting']))
mtext('Density',side=2,line=1.5)
plot.new()
library(MASS)
kde2 <- kde2d(x=faithful[,'waiting'],
              y=faithful[,'eruptions'],
              n=50,lims=c(xlim,ylim))
plot(faithful[,'waiting'],faithful[,'eruptions'],
     col='grey',xlab='waiting',ylab='duration',
     xlim=xlim,ylim=ylim)
contour(kde2,add=TRUE,labcex=0.8,xlim=xlim,ylim=ylim,
        levels=c(0.002,0.005,0.01,0.015,0.02))
mtext('Waiting time [minutes]',side=1,line=1.5)
mtext('Duration [minutes]',side=2,line=1.5)
edens <- density(faithful[,'eruptions'])
plot(edens$y,edens$x,main='',yaxt='n',xlab='',type='l',ylim=ylim)
rug(faithful[,'eruptions'],side=2)
mtext('Density',side=1,line=1.5)
dev.off()

cairo(file='../../figures/ecdfs.pdf',width=8,height=2)
pars(mfrow=c(1,4))
plot(ecdf(pH),main='',verticals=TRUE,pch=NA,
     xlab='pH',ylab='F(pH)')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
plot(ecdf(clasts),main='',verticals=TRUE,pch=NA,
     xlab='clast size (cm)',ylab='F(size)')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
plot(ecdf(por),main='',verticals=TRUE,pch=NA,
     xlab='porosity',ylab='F(porosity)')
legend('topleft',legend='c)',bty='n',cex=1.2,adj=c(2,0))
plot(ecdf(faithful[,'eruptions']),main='',verticals=TRUE,
     pch=NA,xlab='duration [m]',ylab='F(duration)')
legend('topleft',legend='d)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

dat <- pH
xlim <- c(3,7)
mu <- mean(dat)
med <- median(dat)
dens <- density(dat)
mod <- dens$x[which.max(dens$y)]
cairo(file='../../figures/pHlocation.pdf',width=4.5,height=4.5)
pars(mar=rep(0,4))
m <- rbind(c(1,2,3),c(1,4,3),c(1,5,3),c(1,6,3))
layout(m,widths=c(0.1,0.88,0.02),heights=c(0.1,0.4,0.4,0.1))
layout.show(6)
plot.new()
plot.new()
plot.new()
plot(dens,type='l',xaxt='n',main='',xlim=xlim,zero.line=FALSE)
axis(side=3)
mtext(text='pH',side=3,line=1.5)
mtext(text='Density',side=2,line=1.5)
rug(dat,side=1)
lines(rep(mu,2),c(-1,2),lty=1)
lines(rep(med,2),c(-1,2),lty=2)
lines(rep(mod,2),c(-1,2),lty=3)
legend('topleft',legend=c('mean','median','mode'),lty=1:3)
plot(ecdf(dat),main='',verticals=TRUE,pch=NA,xlim=xlim)
mtext(text='pH',side=1,line=1.5)
mtext(text='F(pH)',side=2,line=1.5)
lines(rep(mu,2),c(-1,2),lty=1)
lines(rep(med,2),c(-1,2),lty=2)
lines(rep(mod,2),c(-1,2),lty=3)
lines(xlim,rep(0.5,2),lty=4)
plot.new()
dev.off()

dat <- clasts
xlim <- c(0,12)
mu <- mean(dat)
med <- median(dat)
ldens <- density(log(dat))
dens <- ldens
dens$x <- exp(ldens$x)
dx <- diff(dens$x)
dens$y <- ldens$y/c(dx,dx[511])
mod <- dens$x[which.max(dens$y)]
cairo(file='../../figures/clastslocation.pdf',width=4.5,height=4.5)
pars(mar=rep(0,4))
m <- rbind(c(1,2,3),c(1,4,3),c(1,5,3),c(1,6,3))
layout(m,widths=c(0.1,0.88,0.02),heights=c(0.1,0.4,0.4,0.1))
layout.show(6)
plot.new()
plot.new()
plot.new()
plot(clastdens,type='l',xaxt='n',yaxt='n',main='',xlim=xlim,zero.line=FALSE)
axis(side=2,at=c(0,0.1,0.2,0.3),labels=c(0,0.1,0.2,0.3))
axis(side=3)
mtext(text='clast size [cm]',side=3,line=1.5)
mtext(text='Density',side=2,line=1.5)
rug(dat,side=1)
lines(rep(mu,2),c(-1,30),lty=1)
lines(rep(med,2),c(-1,30),lty=2)
lines(rep(mod,2),c(-1,30),lty=3)
legend('topright',legend=c('mean','median','mode'),lty=1:3)
plot(ecdf(dat),main='',verticals=TRUE,pch=NA,xlim=xlim)
mtext(text='clast size [cm]',side=1,line=1.5)
mtext(text='F(clast size)',side=2,line=1.5)
lines(rep(mu,2),c(-1,2),lty=1)
lines(rep(med,2),c(-1,2),lty=2)
lines(rep(mod,2),c(-1,2),lty=3)
lines(xlim,rep(0.5,2),lty=4)
plot.new()
dev.off()

dat <- por
xlim <- c(0,1)
mu <- mean(dat)
med <- median(dat)
ldens <- density(u)
dens <- ldens
dens$x <- exp(ldens$x)/(exp(ldens$x)+1)
dx <- diff(dens$x)
dens$y <- ldens$y/c(dx,dx[511])
mod <- dens$x[which.max(dens$y)]
cairo(file='../../figures/porositylocation.pdf',width=4.5,height=4.5)
pars(mar=rep(0,4))
m <- rbind(c(1,2,3),c(1,4,3),c(1,5,3),c(1,6,3))
layout(m,widths=c(0.1,0.88,0.02),heights=c(0.1,0.4,0.4,0.1))
layout.show(6)
plot.new()
plot.new()
plot.new()
plot(pordens,type='l',xaxt='n',main='',xlim=xlim,zero.line=FALSE)
axis(side=3)
mtext(text='Porosity',side=3,line=1.5)
mtext(text='Density',side=2,line=1.5)
mtext(text='a)',side=3,at=-0.11,line=1.5,cex=1.2)
legend(0.1,114,legend=c('mean','median','mode'),lty=1:3)
rug(dat,side=1)
lines(rep(mu,2),c(-1,120),lty=1)
lines(rep(med,2),c(-1,120),lty=2)
lines(rep(mod,2),c(-1,120),lty=3)
plot(ecdf(dat),main='',verticals=TRUE,pch=NA,xlim=xlim)
mtext(text='Porosity',side=1,line=1.5)
mtext(text='F(porosity)',side=2,line=1.5)
lines(rep(mu,2),c(-1,2),lty=1)
lines(rep(med,2),c(-1,2),lty=2)
lines(rep(mod,2),c(-1,2),lty=3)
lines(xlim,rep(0.5,2),lty=4)
plot.new()
dev.off()

dat <- faithful[,'eruptions']
xlim <- c(0.8,6)
mu <- mean(dat)
med <- median(dat)
dens <- density(dat)
mod <- dens$x[which.max(dens$y)]
cairo(file='../../figures/eruptionslocation.pdf',width=4.5,height=4.5)
pars(mar=rep(0,4))
m <- rbind(c(1,2,3),c(1,4,3),c(1,5,3),c(1,6,3))
layout(m,widths=c(0.1,0.88,0.02),heights=c(0.1,0.4,0.4,0.1))
plot.new()
plot.new()
plot.new()
plot(dens,type='l',xaxt='n',main='',xlim=xlim,zero.line=FALSE)
axis(side=3)
mtext(text='Duration [min]',side=3,line=1.5)
mtext(text='Density',side=2,line=1.5)
mtext(text='b)',side=3,adj=-0.09,line=1.5,cex=1.2)
rug(dat,side=1)
lines(rep(mu,2),c(-1,2),lty=1)
lines(rep(med,2),c(-1,2),lty=2)
lines(rep(mod,2),c(-1,2),lty=3)
legend('topleft',legend=c('mean','median','mode'),lty=1:3)
plot(ecdf(dat),main='',verticals=TRUE,pch=NA,xlim=xlim)
mtext(text='Duration [min]',side=1,line=1.5)
mtext(text='F(duration)',side=2,line=1.5)
lines(rep(mu,2),c(-1,2),lty=1)
lines(rep(med,2),c(-1,2),lty=2)
lines(rep(mod,2),c(-1,2),lty=3)
lines(xlim,rep(0.5,2),lty=4)
plot.new()
dev.off()

cairo(file='../../figures/IQRpH.pdf',width=4.5,height=4.5)
pars(mar=c(2.5,2.3,0.5,0.4))
plot(ecdf(pH),main='',col.01line=NULL,verticals=TRUE,
     pch=NA,xlab='pH',ylab='F(pH)')
quants <- quantile(pH,c(0.25,0.5,0.75))
lines(c(2,quants['25%'],quants['25%']),c(0.25,0.25,-0.1),lty=2)
lines(c(2,quants['50%'],quants['50%']),c(0.5,0.5,-0.1),lty=4)
lines(c(2,quants['75%'],quants['75%']),c(0.75,0.75,-0.1),lty=3)
legend('topleft',lty=c(2,4,3),bg='white',
       legend=c('25 percentile','median','75 percentile'))
dev.off()

skewness <- function(x){
    mean((x-mean(x))^3)/sd(x)^3
}
covid <- read.csv('covid.csv',header=TRUE,check.names=FALSE)
age <- (covid[,'to']+covid[,'from'] )/2
meanage <- sum(age*covid[,'death rate'])/sum(covid[,'death rate'])
stdage <- sqrt(sum(covid[,'death rate']*(age-meanage)^2)/sum(covid[,'death rate']))
skewage <- sum(covid[,'death rate']*(age-meanage)^3)/
    (sum(covid[,'death rate'])*stdage^3)
cairo(file='../../figures/skewness.pdf',width=7,height=2)
pars(mar=c(2.5,2.3,1,0.2),mfrow=c(1,3))
plot(density(pH),main='',xlab='pH',bty='n',zero.line=FALSE)
mtext(paste0('skewness=',signif(skewness(pH),2)),side=3,line=0,cex=0.8)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
ldens <- density(log(clasts))
dens <- ldens
dens$x <- exp(ldens$x)
dx <- diff(dens$x)
dens$y <- ldens$y/c(dx,dx[511])
plot(clastdens$x,clastdens$y,main='',xlab='clast size [cm]',
     ylab='Density',type='l',xlim=c(0,12),bty='n')
mtext(paste0('skewness=',signif(skewness(clasts),2)),side=3,line=0,cex=0.8)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(0,0))
barplot(height=covid[,'death rate'],
        width=covid[,'to']-covid[,'from']+1,
        space=0,col=NA,xlim=c(0,100),ylim=c(0,1200),
        xlab='age of death [years]',ylab='deaths/100k people')
mtext(paste0('skewness=',signif(skewage,2)),side=3,line=0,cex=0.8)
axis(age,side=1,at=c(0,25,50,75,100),labels=c(0,25,50,75,100))
legend('topleft',legend='c)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

cairo(file='../../figures/boxplot.pdf',width=4,height=3)
pars(mar=rep(0,4))
m <- rbind(c(1,2),c(1,3),c(4,4))
layout(m,widths=c(0.1,0.9),heights=c(0.65,0.2,0.15))
plot.new()
xlim <- c(0,12)
iin <- which(dens$x<xlim[2])
plot(clastdens$x[iin],clastdens$y[iin],type='l',
     main='',xlim=xlim,xaxt='n',bty='n',ylab='')
mtext('Density',side=2,line=1.5)
rug(clasts,ticksize=0.05)
boxplot(clasts,horizontal=TRUE,ylim=xlim,
        frame.plot=FALSE,xlab='',width=1)
mtext('clast size [cm]',side=1,line=1.5)
dev.off()

nn <- 5
H0 <- 2/3

cairo(file='../../figures/goldbarplot.pdf',width=3,height=3)
pars()
prob <- dbinom(0:nn,nn,H0)
names(prob) <- 0:nn
barplot(prob,col=NA,xlab='k = # gold discoveries',ylab='P(k)',space=0)
dev.off()

cairo(file='../../figures/goldCDF.pdf',width=3,height=3)
pars()
X <- -1:(nn+1)
P <- pbinom(X,nn,H0)
plot(X,P,type='s',xlab='x = # gold discoveries',
     ylab='F(x)',bty='n',xaxt='n')
axis(side=1,at=c(0:nn))
dev.off()

binomhist <- function(nn,kk,H0,Ha=H0,nsides=1,rej.col='black',
                      na.col=NA, showax=TRUE,plotk=TRUE,
                      xlab='k = # gold discoveries',ylab='P(k)',...){
    alpha <- 0.05
    prob <- dbinom(0:nn,nn,Ha)
    names(prob) <- 0:nn
    if (nsides==-1){
        lrej <- 0
        urej <- qbinom(0.95,nn,H0)
    } else if (nsides==1){
        lrej <- qbinom(0.05,nn,H0)
        urej <- 0
    } else {
        lrej <- qbinom(0.025,nn,H0)
        urej <- nn-qbinom(0.975,nn,H0)
    }
    nacc <- nn+1-lrej-urej
    if (showax){
        barplot(prob,col=c(rep(rej.col,lrej),rep(na.col,nacc),rep(rej.col,urej)),
                xlab=xlab,ylab=ylab,space=0,...)
    } else {
        barplot(prob,col=c(rep(rej.col,lrej),rep(na.col,nacc),rep(rej.col,urej)),
                space=0,xlab='',ylab='',xaxt='n',yaxt='n',...)
    }
    if (plotk) lines(rep(kk,2)+0.5,c(0,1),lty=2)
}

binomcdf <- function(nn,kk,H0,Ha=H0,nsides=1,showax=TRUE,add=FALSE,
                     plotk=TRUE,plota=TRUE,plotp=TRUE,dist='binomial',
                     xlab='x = # gold discoveries',ylab='F(x)',...){
    X <- -1:(nn+1)
    P <- pbinom(X,nn,Ha)
    if (add){
        lines(X,P,type='s',...)
    } else if (showax){
        plot(X,P,type='s',xlab=xlab,ylab=ylab,bty='n',xaxt='n',...)
        axis(side=1,at=c(0:nn))
    } else {
        plot(X,P,type='s',xlab='',ylab='',bty='n',xaxt='n',yaxt='n',...)
    }
    if (!add){
        if (plota){
            if (nsides%in%c(-1,1)){
                if (nsides==1){
                    lines(c(-1,nn+1),rep(0.05,2),lty=3)
                } else {
                    lines(c(-1,nn+1),rep(0.95,2),lty=3)
                }
                if (plotp){
                    pval <- pbinom(kk,nn,H0)
                    lines(c(-1,nn+1),rep(pval,2),lty=2)
                } else {
                    foo <- 1
                }
            } else {
                lines(c(-1,nn+1),rep(0.025,2),lty=3)
                lines(c(-1,nn+1),rep(0.975,2),lty=3)
                if (plotp){
                    pval <- pbinom(kk,nn,H0)
                    lines(c(-1,nn+1),rep(pval,2),lty=2)
                } else {
                    foo <- 1
                }
            }
        }
        if (nsides%in%c(-1,1)){
            if (nsides==1){
                lrej <- qbinom(0.05,nn,H0)
            } else {
                lrej <- qbinom(0.95,nn,H0)
            }
            urej <- 0
            lines(rep(lrej,2),c(0,1),lty=3)
        } else if (nsides==2){
            lrej <- qbinom(0.025,nn,H0)
            urej <- nn-qbinom(0.975,nn,H0)
            lines(rep(lrej,2),c(0,1),lty=3)
            lines(rep(nn-urej,2),c(0,1),lty=3)
        }
        if (plotk) lines(rep(kk,2),c(0,1),lty=2)
    }
}

cairo(file='../../figures/1sidedbinomialrejection5.pdf',width=6,height=3)
pars(mfrow=c(1,2))
binomhist(nn=5,kk=2,H0=2/3,nsides=1)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=5,kk=2,H0=2/3,nsides=1)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

cairo(file='../../figures/2sidedbinomialrejection5.pdf',width=6,height=3)
pars(mfrow=c(1,2))
binomhist(nn=5,kk=2,H0=2/3,nsides=2)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=5,kk=2,H0=2/3,nsides=2)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,1))
dev.off()

cairo(file='../../figures/1sidedbinomialrejection15.pdf',width=6,height=3)
pars(mfrow=c(1,2))
binomhist(nn=15,kk=6,H0=2/3,nsides=1)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=15,kk=6,H0=2/3,nsides=1)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

cairo(file='../../figures/2sidedbinomialrejection15.pdf',width=6,height=3)
pars(mfrow=c(1,2))
binomhist(nn=15,kk=6,H0=2/3,nsides=2)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=15,kk=6,H0=2/3,nsides=2)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,1))
dev.off()

cairo(file='../../figures/binomialrejection30.pdf',width=5,height=3.5)
pars(mar=rep(0,4))
m <- rbind(c(1,2,3),c(1,4,5),c(1,6,6))
layout(m,widths=c(0.1,0.45,0.45),heights=c(0.4,0.4,0.2))
plot.new()
binomhist(nn=30,kk=12,H0=2/3,nsides=1,showax=FALSE,xlim=c(-0.5,31.5))
axis(side=2); mtext('P(k)',side=2,line=2,cex=0.8)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomhist(nn=30,kk=12,H0=2/3,nsides=2,showax=FALSE,xlim=c(-0.5,31.5))
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=30,kk=12,H0=2/3,nsides=1,showax=FALSE,xlim=c(-1,31))
axis(side=1); mtext('# gold discoveries',cex=0.8,side=1,line=2)
axis(side=2); mtext('F(x)',side=2,cex=0.8,line=2)
legend('topleft',legend='c)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=30,kk=12,H0=2/3,nsides=2,showax=FALSE,xlim=c(-1,31))
axis(side=1); mtext('# gold discoveries',cex=0.8,side=1,line=2)
legend('topleft',legend='d)',bty='n',cex=1.2,adj=c(2,1))
plot.new()
dev.off()

cairo(file='../../figures/binompower1.pdf',width=8,height=3.8)
pars(mar=rep(0,4))
m <- rbind(c(1,10,10,10),c(1,2,3,4),c(1,8,8,8),c(1,5,6,7),c(1,9,9,9))
layout(m,widths=c(0.05,0.33,0.33,0.33),
       heights=c(0.02,0.42,0.04,0.42,0.1))
plot.new()
binomhist(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),ylim=c(0,1),border='white',
          rej.col='gray70',na.col='gray70',plotk=FALSE)
binomhist(nn=5,kk=2,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),na.col='white',add=TRUE,plotk=FALSE)
axis(side=2); mtext('P(k)',side=2,cex=0.8,line=1.5)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomhist(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),ylim=c(0,1),border='white',
          rej.col='gray70',na.col='gray70',plotk=FALSE)
binomhist(nn=5,kk=2,H0=2/3,Ha=1/5,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),na.col='white',add=TRUE,plotk=FALSE)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
binomhist(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),ylim=c(0,1),border='white',
          rej.col='gray70',na.col='gray70',plotk=FALSE)
binomhist(nn=5,kk=2,H0=2/3,Ha=0,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),na.col='white',add=TRUE,plotk=FALSE)
legend('topleft',legend='c)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
         xlim=c(-1,6),col='gray60',plotp=FALSE,plotk=FALSE)
binomcdf(nn=5,kk=2,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
         xlim=c(-1,6),add=TRUE,plotp=FALSE,plotk=FALSE)
b <- pbinom(qbinom(0.05,5,2/3)-1,5,2/5)
lines(c(-1,6),rep(b,2),lty=2)
axis(side=1,at=0:5); mtext('# gold discoveries',side=1,cex=0.8,line=1.5)
axis(side=2)
mtext('P(k)',side=2,cex=0.8,line=1.5)
legend('topleft',legend='d)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
         xlim=c(-1,6),col='gray60',plotp=FALSE,plotk=FALSE)
binomcdf(nn=5,kk=2,H0=2/3,Ha=1/5,nsides=1,showax=FALSE,
         xlim=c(-1,6),add=TRUE,plotp=FALSE,plotk=FALSE)
b <- pbinom(qbinom(0.05,5,2/3)-1,5,1/5)
lines(c(-1,6),rep(b,2),lty=2)
axis(side=1,at=0:5)
mtext('# gold discoveries',side=1,cex=0.8,line=1.5)
legend('topleft',legend='e)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
         xlim=c(-1,6),col='gray60',plotp=FALSE,plotk=FALSE)
binomcdf(nn=5,kk=2,H0=2/3,Ha=0,nsides=1,plotk=FALSE,
         showax=FALSE,xlim=c(-1,6),plotp=FALSE,add=TRUE)
b <- pbinom(qbinom(0.05,5,2/3)-1,5,0)
lines(c(-1,16),rep(b,2),lty=2)
axis(side=1,at=0:5)
mtext('# gold discoveries',side=1,cex=0.8,line=1.5)
legend('topleft',legend='f)',bty='n',cex=1.2,adj=c(2,1))
plot.new()
plot.new()
dev.off()

cairo(file='../../figures/binompower2.pdf',width=8,height=3.8)
pars(mar=rep(0,4))
m <- rbind(c(1,2,3,4),c(1,8,8,8),c(1,5,6,7),c(1,9,9,9))
layout(m,widths=c(0.05,0.33,0.33,0.33),
       heights=c(0.43,0.04,0.43,0.1))
plot.new()
binomhist(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),border='white',rej.col='gray70',
          na.col='gray70',ylim=c(0,0.35),plotk=FALSE)
binomhist(nn=5,kk=2,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
          xlim=c(-0.5,6.5),na.col='white',add=TRUE,plotk=FALSE)
axis(side=2); mtext('P(k)',side=2,cex=0.8,line=1.5)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
binomhist(nn=15,kk=6,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
          xlim=c(-0.5,16.5),border='white',
          rej.col='gray70',na.col='gray70',plotk=FALSE)
binomhist(nn=15,kk=6,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
          xlim=c(-0.5,16.5),na.col='white',add=TRUE,plotk=FALSE)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
binomhist(nn=30,kk=12,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
          xlim=c(-0.5,31.5),border='white',
          rej.col='gray70',na.col='gray70',plotk=FALSE)
binomhist(nn=30,kk=12,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
          xlim=c(-0.5,31.5),na.col='white',add=TRUE,plotk=FALSE)
legend('topleft',legend='c)',bty='n',cex=1.2,adj=c(2,0))
binomcdf(nn=5,kk=2,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
         xlim=c(-1,6),col='gray60',plotp=FALSE,plotk=FALSE)
binomcdf(nn=5,kk=2,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
         xlim=c(-1,6),add=TRUE,plotp=FALSE,plotk=FALSE)
b <- pbinom(qbinom(0.05,5,2/3)-1,5,2/5)
lines(c(-1,6),rep(b,2),lty=2)
axis(side=1,at=0:5); mtext('# gold discoveries',side=1,cex=0.8,line=1.5)
axis(side=2)
mtext('P(k)',side=2,cex=0.8,line=1.5)
legend('topleft',legend='d)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=15,kk=6,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
         xlim=c(-1,16),col='gray60',plotp=FALSE,plotk=FALSE)
binomcdf(nn=15,kk=6,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
         xlim=c(-1,16),add=TRUE,plotp=FALSE,plotk=FALSE)
b <- pbinom(qbinom(0.05,15,2/3)-1,15,2/5)
lines(c(-1,16),rep(b,2),lty=2)
axis(side=1,at=seq(from=0,to=15,by=5))
mtext('# gold discoveries',side=1,cex=0.8,line=1.5)
legend('topleft',legend='e)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=30,kk=12,H0=2/3,Ha=2/3,nsides=1,showax=FALSE,
         xlim=c(-1,31),col='gray60',plotp=FALSE,plotk=FALSE)
binomcdf(nn=30,kk=12,H0=2/3,Ha=2/5,nsides=1,showax=FALSE,
         xlim=c(-1,31),add=TRUE,plotp=FALSE,plotk=FALSE)
b <- pbinom(qbinom(0.05,30,2/3)-1,30,2/5)
lines(c(-1,31),rep(b,2),lty=2)
axis(side=1,at=seq(from=0,to=30,by=5))
mtext('# gold discoveries',side=1,cex=0.8,line=1.5)
legend('topleft',legend='f)',bty='n',cex=1.2,adj=c(2,1))
plot.new()
plot.new()
dev.off()

cairo(file='../../figures/binompvsn.pdf',width=3.7,height=3)
pars()
k <- c(2*c(2,5,10)-1,2*c(2,5,10)+1)
n <- rep(3*c(2,5,10),2)
lty <- rep(c(1,2,4),2)
col <- c(rep('black',3),rep('gray60',3))
maxn <- 1000
plot(c(0,maxn),c(0,.4),type='n',xlab='n',ylab='p-value',bty='n')
leg <- list()
for (i in 1:length(k)){
    leg <- c(leg,paste0('k/n=',k[i],'/',n[i],' (=',signif(k[i]/n[i],3),')'))
    multiplier <- 1:floor(maxn/n[i])
    pval <- pbinom(k[i]*multiplier,n[i]*multiplier,2/3)
    if (any(pval>0.5)) pval <- 1-pval
    lines(n[i]*multiplier,pval,col=col[i],lty=lty[i])
}
lines(c(0,maxn),rep(0.05,2),lty=3)
legend('topright',legend=leg,lty=lty,col=col,cex=0.8)
dev.off()

cairo(file='../../figures/binomcik2n5.pdf',width=7,height=3.5)
n <- 5
k <- 2
pars(mfrow=c(1,2))
ci <- binom.test(k,n)$conf.int
binomhist(nn=n,kk=k,H0=ci[1],Ha=ci[1],nsides=2,
          border='white',na.col='gray60')
binomhist(nn=n,kk=k,H0=ci[2],Ha=ci[2],nsides=2,add=TRUE)
legend('topright',legend='a)',bty='n',cex=1.2,adj=c(0,1))
binomcdf(nn=n,kk=k,H0=ci[1],Ha=ci[1],nsides=0,
         plota=TRUE,plotp=FALSE,col='gray60')
ly <- pbinom(0,n,ci[1])
text(0,ly,labels=paste0('p=',signif(ci[1],2)),pos=3,col='gray60')
binomcdf(nn=n,kk=k,H0=ci[2],Ha=ci[2],nsides=0,add=TRUE)
ly <- pbinom(4,n,ci[2])
text(5,ly,labels=paste0('p=',signif(ci[2],2)),pos=1)
legend('topright',legend='b)',bty='n',cex=1.2,adj=c(0,1))
dev.off()

cairo(file='../../figures/binomcik4n5.pdf',width=7,height=3.5)
n <- 5
k <- 4
pars(mfrow=c(1,2))
ci <- binom.test(k,n)$conf.int
binomhist(nn=n,kk=k,H0=ci[1],Ha=ci[1],nsides=2,
          ylim=c(0,1),border='white',na.col='gray60')
binomhist(nn=n,kk=k,H0=ci[2],Ha=ci[2],nsides=2,add=TRUE)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=n,kk=k,H0=ci[1],Ha=ci[1],nsides=0,
         plota=TRUE,plotp=FALSE,col='gray60')
ly <- pbinom(0,n,ci[1])
text(0,ly,labels=paste0('p=',signif(ci[1],2)),pos=3,col='gray60')
binomcdf(nn=n,kk=k,H0=ci[2],Ha=ci[2],nsides=0,add=TRUE)
ly <- pbinom(4,n,ci[2])
text(3,ly,labels=paste0('p=',signif(ci[2],2)),pos=1,offset=-0.75)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,1))
dev.off()

cairo(file='../../figures/binomcik12n30.pdf',width=7,height=3.5)
n <- 30
k <- 12
pars(mfrow=c(1,2))
ci <- binom.test(k,n)$conf.int
binomhist(nn=n,kk=k,H0=ci[1],Ha=ci[1],nsides=2,showax=FALSE,
          border='white',na.col='gray60')
binomhist(nn=n,kk=k,H0=ci[2],Ha=ci[2],
          showax=FALSE,nsides=2,add=TRUE)
axis(side=1,at=seq(from=0,to=30,by=5))
axis(side=2)
mtext('# gold discoveries',side=1,line=1.5)
mtext('P(k)',side=2,line=1.5)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1))
binomcdf(nn=n,kk=k,H0=ci[1],Ha=ci[1],nsides=0,
         showax=FALSE,plota=TRUE,plotp=FALSE,col='gray60')
ly <- pbinom(6,n,ci[1])
text(6,ly,labels=paste0('p=',signif(ci[1],2)),pos=2,col='gray60')
binomcdf(nn=n,kk=k,H0=ci[2],Ha=ci[2],showax=FALSE,nsides=0,add=TRUE)
ly <- pbinom(17,n,ci[2])
text(17,ly,labels=paste0('p=',signif(ci[2],2)),pos=4)
axis(side=1,at=seq(from=0,to=30,by=5))
axis(side=2)
mtext('# gold discoveries',side=1,line=1.5)
mtext('F(x)',side=2,line=1.5)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,1))
dev.off()

cairo(file='../../figures/binomcivsn.pdf',width=7,height=3)
pars(mfrow=c(1,2))
k <- 2
n <- 3
multiplier <- c(1,seq(from=10,to=100,by=10))
maxn <- max(multiplier)*n
plot(c(0,maxn),c(0,1),type='n',xlab='n',ylab='p',bty='n')
for (m in multiplier){
    ci <- binom.test(k*m,n*m,p=k/n)$conf.int
    arrows(n*m,ci[1],n*m,ci[2],code=3,angle=90,length=0.05)
}
lines(c(0,maxn),rep(k/n,2),lty=3)
legend('topright',legend='a)',bty='n',cex=1.2,adj=c(0,0))
k <- 1
n <- 5
multiplier <- c(1,seq(from=5,to=60,by=5))
maxn <- max(multiplier)*n
plot(c(0,maxn),c(0,1),type='n',xlab='n',ylab='p',bty='n')
for (m in multiplier){
    ci <- binom.test(k*m,n*m,p=k/n)$conf.int
    arrows(n*m,ci[1],n*m,ci[2],code=3,angle=90,length=0.05)
}
legend('topright',legend='b)',bty='n',cex=1.2,adj=c(0,0))
lines(c(0,maxn),rep(k/n,2),lty=3)
dev.off()

cairo(file='../../figures/declusteredquakes.pdf',width=8,height=1.8)
pars()
dat <- read.csv('declusteredquakes.csv',header=TRUE)
from <- 1917
to <- 2016
minmag <- 5.0
bigrecent <- dat$Year[which(dat$Year>=from & dat$Year<=to & dat$Mwe>=minmag)]
tab <- table(bigrecent)
years <- from:to
nquakes <- rep(0,length(years))
names(nquakes) <- years
nquakes[names(tab)] <- tab
b <- barplot(nquakes,space=1,col='black',
             xlab='year',ylab='# earthquakes',xaxt='n')
ticks <- unique(c(seq.int(from=from,to=to,by=5),to))
axis(side=1,at=b[1]+(ticks-from)*(b[length(b)]-b[1])/(to-from),labels=ticks)
dev.off()

cairo(file='../../figures/declusteredquakesperyear.pdf',width=4,height=3)
pars()
hist(nquakes,breaks=seq(from=-0.5,to=12.5,by=1),main='',
     xlab='# earthquakes/year')
dev.off()

cairo(file='../../figures/zircons.pdf',width=15,height=5)
pars(mar=c(1,0,0,0))
nx <- 12
ny <- 4
nn <- 5000/(nx*ny)
pd <- 0.03
seeds <- (1:(4*nx*ny+2))+22
set.seed(seeds[1])
nd <- matrix(rpois(nx*ny,pd*nn),nx,ny)
set.seed(seeds[2])
no <- matrix(rpois(nx*ny,(1-pd)*nn),nx,ny)
boxx <- c(0,0,nx,nx,0)
boxy <- c(0,ny,ny,0,0)
hlinex <- rbind(rep(0,ny+1),rep(nx,ny+1))
hliney <- rbind(0:ny,0:ny)
vlinex <- rbind(0:nx,0:nx)
vliney <- rbind(rep(0,nx+1),rep(ny,nx+1))
plot(boxx,boxy,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
polygon(boxx,boxy)
seednum <- -1
for (x in 1:nx){
    for (y in 1:ny){
        seednum <- seednum+4
        set.seed(seeds[seednum])
        xo <- x + runif(no[x,y]) - 1
        set.seed(seeds[seednum+1])
        yo <- y + runif(no[x,y]) - 1
        set.seed(seeds[seednum+2])
        xd <- x + runif(nd[x,y]) - 1
        set.seed(seeds[seednum+3])
        yd <- y + runif(nd[x,y]) - 1
        points(xo,yo,pch=21,bg='white',col='gray60')
        points(xd,yd,pch=22,fg='black',bg='black')
    }
}
matlines(hlinex,hliney,col='gray30',lty=2)
matlines(vlinex,vliney,col='gray30',lty=2)
p <- par(xpd=TRUE)
text(x=seq(from=0.5,to=nx-0.5,by=1),y=-0.2,labels=1:nx,cex=1.5,col='gray30')
text(x=-0.2,y=seq(from=0.5,to=ny-0.5,by=1),labels=1:ny,cex=1.5,col='gray30')
par(p)
dev.off()

cairo(file='../../figures/zirconcounts.pdf',width=12,height=4)
pars(mar=c(1,0,0,0))
plot(boxx,boxy,type='n',bty='n',xaxt='n',yaxt='n',xlab='',ylab='')
polygon(boxx,boxy)
matlines(hlinex,hliney,col='gray30',lty=2)
matlines(vlinex,vliney,col='gray30',lty=2)
p <- par(xpd=TRUE)
text(expand.grid(1:nx,1:ny)-0.5,labels=nd)
text(x=seq(from=0.5,to=nx-0.5,by=1),y=-0.2,labels=1:nx,col='gray30')
text(x=-0.2,y=seq(from=0.5,to=ny-0.5,by=1),labels=1:ny,col='gray30')
par(p)
dev.off()

cairo(file='../../figures/zirconhist.pdf',width=3,height=2.5)
pars()
nzirc <- as.vector(nd)
h <- table(nzirc)
barplot(h,col='white',xlab='# zircons/graticule',ylab='# graticules')
dev.off()

# nn is a vector of values to be evaluated
poishist <- function(nn,kk,H0,Ha=H0,nsides=1,
                     rej.col='black',na.col=NA,
                     showax=TRUE,plotk=TRUE,plotl=TRUE,
                     xlab='k',ylab='P(k)',...){
    alpha <- 0.05
    prob <- dpois(nn,Ha)
    names(prob) <- nn
    if (nsides==-1){
        lrej <- 0
        urej <- length(nn)-qpois(0.95,H0)+1
    } else if (nsides==1){
        lrej <- qpois(0.05,H0)
        urej <- 0
    } else {
        lrej <- qpois(0.025,H0)
        urej <- length(nn)-qpois(0.975,H0)
    }
    nacc <- length(nn)+1-lrej-urej
    if (showax){
        barplot(prob,col=c(rep(rej.col,lrej),rep(na.col,nacc),rep(rej.col,urej)),
                xlab=xlab,ylab=ylab,space=0,...)
    } else {
        barplot(prob,col=c(rep(rej.col,lrej),rep(na.col,nacc),rep(rej.col,urej)),
                space=0,xlab='',ylab='',xaxt='n',yaxt='n',...)
    }
    if (plotk) lines(rep(kk,2)+0.5,c(0,1),lty=2)
    if (plotl) lines(rep(H0,2)+0.5,c(0,1),lty=2)
}

poiscdf <- function(nn,kk,H0,Ha=H0,nsides=1,showax=TRUE,add=FALSE,
                    plotk=TRUE,plota=TRUE,plotp=TRUE,plotl=TRUE,
                    xlab='k',ylab='F(k)',...){
    X <- nn
    P <- ppois(X,Ha)
    if (add){
        lines(X,P,type='s',...)
    } else if (showax){
        plot(X,P,type='s',xlab=xlab,ylab=ylab,bty='n',xaxt='n',...)
        axis(side=1,at=nn)
    } else {
        plot(X,P,type='s',xlab='',ylab='',bty='n',xaxt='n',yaxt='n',...)
    }
    if (!add){
        if (plota){
            if (nsides%in%c(-1,1)){
                if (nsides==1){
                    lines(range(nn),rep(0.05,2),lty=3)
                } else {
                    lines(range(nn),rep(0.95,2),lty=3)
                }
                if (plotp){
                    pval <- ppois(kk,H0)
                    lines(range(nn),rep(pval,2),lty=2)
                }
            } else {
                lines(range(nn),rep(0.025,2),lty=3)
                lines(range(nn),rep(0.975,2),lty=3)
                if (plotp){
                    pval <- ppois(kk,H0)
                    lines(range(nn),rep(pval,2),lty=2)
                }
            }
        }
        if (nsides%in%c(-1,1)){
            if (nsides==1){
                lrej <- qpois(0.05,H0)
                urej <- 0
                lines(rep(lrej,2),c(0,1),lty=3)
            } else {
                lrej <- 0
                urej <- qpois(0.95,H0)
                lines(rep(urej,2),c(0,1),lty=3)
            }
        } else if (nsides==2){
            lrej <- qpois(0.025,H0)
            urej <- max(nn)-qpois(0.975,H0)
            lines(rep(lrej,2),c(0,1),lty=3)
            lines(rep(max(nn)-urej,2),c(0,1),lty=3)
        }
        if (plotk) lines(rep(kk,2),c(0,1),lty=2)
        if (plotl) lines(rep(H0,2),c(0,1),lty=2)
    }
}

cairo(file='../../figures/increasinglambda.pdf',width=8,height=3)
lambda <- c(1,2,5,10)
labels <- c('a)','b)','c)','d)')
pars(mar=rep(0,4),mfrow=c(2,length(lambda)))
m <- rbind(c(1,rep(10,7)),
           c(1,2,13,3,14,4,15,5),
           c(1,rep(11,7)),
           c(1,6,16,7,17,8,18,9),
           c(1,rep(12,7)))
layout(m,widths=c(0.05,0.22,0.02,0.22,0.02,0.22,0.02,0.22),
       heights=c(0.05,0.39,0.05,0.39,0.12))
plot.new()
for (i in 1:length(lambda)){
    M <- qpois(0.99999,lambda[i])
    poishist(0:M,H0=lambda[i],plotk=FALSE,rej.col='white',
             showax=TRUE,xlim=c(-1,M)+0.5,xlab='',xaxt='n',yaxt='n')
    axis(side=2,pos=-(M/50))
    if (i==1) mtext(side=2,line=1.5,text='P(k)',cex=0.9)
    legend('topright',legend=labels[i],bty='n',cex=1.5,adj=c(1,0))
}
labels <- c('e)','f)','g)','h)')
for (i in 1:length(lambda)){
    M <- qpois(0.99999,lambda[i])
    poiscdf(0:M,H0=lambda[i],plotk=FALSE,plotp=FALSE,
            plota=FALSE,showax=FALSE,nsides=0,xlim=c(-1,M))
    if (i==1) {
        axis(side=2,pos=-0.5-(M/50))
        mtext(side=2,line=1.5,text='F(x)',cex=0.9)
    }
    axis(side=1,at=pretty(0:M))
    mtext(side=1,line=1.5,text='x,k',cex=0.9)
    legend('topleft',legend=labels[i],bty='n',cex=1.5,adj=c(1,0))
}
plot.new()
dev.off()

if (FALSE){
    lambda <- 5
    nn <- c(10,20,50,100,200,500,1000,2000,5000,10000)
    p <- lambda/nn
    quant <- qpois(0.1,lambda)
    pb <- pbinom(quant,nn,p)
    pp <- ppois(quant,lambda)
    signif(rbind(nn,p,pb),3)
}

cairo(file='../../figures/poishyp.pdf',width=6,height=3)
pars(mfrow=c(1,2))
nn <- 0:10
l <- 3.5
poishist(nn,H0=l,kk=9,plotk=TRUE,rej.col='black',showax=TRUE,nsides=-1)
poiscdf(nn,H0=l,kk=9,plotk=TRUE,plotp=TRUE,plotl=FALSE,
        plota=TRUE,showax=TRUE,nsides=-1)
dev.off()

cairo(file='../../figures/poisci.pdf',width=6,height=2.7)
pars(mfrow=c(1,2))
nn <- 0:25
kk <- 5
ci <- poisson.test(kk)$conf.int
poishist(nn,H0=ci[1],kk=kk,na.col='gray60',
         plotk=TRUE,plotl=FALSE,
         rej.col='black',border='white',showax=TRUE,nsides=2)
poishist(nn,H0=ci[2],kk=kk,rej.col='black',
         plotk=FALSE,plotl=FALSE,showax=FALSE,nsides=2,add=TRUE)
poiscdf(nn,H0=ci[1],kk=kk,col='gray60',
        plotk=TRUE,plota=TRUE,plotp=FALSE,plotl=FALSE,
        showax=TRUE,nsides=0,ylim=c(0,1))
text(ci[1],0.15,col='gray60',pos=4,srt=90,
     labels=bquote(lambda == .(signif(ci[1],3))))
poiscdf(nn,H0=ci[2],kk=kk,add=TRUE,nsides=0,
        plotk=FALSE,plota=FALSE,plotp=FALSE,plotl=FALSE)
text(ci[2],0.15,col='black',pos=4,srt=90,
     labels=bquote(lambda == .(signif(ci[2],3))))
dev.off()

cairo(file='../../figures/poiserrbars.pdf',width=10,height=2.5)
pars(mar=c(2.5,2.6,0.0,0.0))
plot(range(years),c(0,21),type='n',bty='n',
     xlab='year',ylab=expression('# earthquakes/year ('*lambda*')'),xaxt='n')
axis(side=1,at=c(1917,1925,1950,1975,2000,2016))
lambda <- mean(nquakes)
for (i in 1:length(nquakes)){
    pt <- poisson.test(nquakes[i],r=lambda)
    ci <- pt$conf.int
    if (ci[1]<lambda & ci[2]>lambda) col <- 'gray60'
    else col <- 'black'
    points(years[i],nquakes[i],pch=20,col=col,cex=0.8)
    arrows(years[i],ci[1],years[i],ci[2],
           code=3,angle=90,length=0.02,col=col)
}
lines(range(years),rep(mean(nquakes),2),lty=3)
dev.off()

cairo(file='../../figures/quakerug.pdf',width=8,height=1)
dat <- read.csv('declusteredquakes.csv',header=TRUE)
pars(mar=c(2.5,0,0,0))
from <- 1917
to <- 2016
minmag <- 5.0
bigrecent <- (dat$Year>=from & dat$Mwe>=minmag)
allzero <- (dat$Hour==0 & dat$Minute==0 & dat$Second==0)
good <- bigrecent & !allzero
eqtimes <- dat$Hour[good] + dat$Minute[good]/60 + dat$Second[good]/3600
plot(eqtimes,rep(0,length(eqtimes)),type='p',bty='n',pch='|',
     yaxt='n',ylab='',xaxt='n',xlab='time (hr)')
axis(side=1,at=seq(from=0,to=24,by=2))
dev.off()

CLT <- function(n,dat,xlab,plot=TRUE){
    m <- 500
    s <- rep(0,m)
    for (i in 1:m){
        samp <- sample(dat,n)
        s[i] <- sum(samp)
    }
    dens <- density(s)
    if (plot){
        plot(dens,main='',zero.line=FALSE,xlab=xlab,bty='n')
        rug(s)
    }
    invisible(s)
}

cairo(file='../../figures/CLTfaithful1.pdf',width=3,height=2)
pars()
dens <- density(faithful[,'eruptions'])
plot(dens,main='',zero.line=FALSE,xlab='duration (min)',bty='n')
rug(faithful[,'eruptions'])
dev.off()

set.seed(7)

cairo(file='../../figures/CLTfaithful2.pdf',width=3,height=2)
pars()
CLT(2,dat=faithful[,'eruptions'],xlab='total duration (min)')
dev.off()

cairo(file='../../figures/CLTfaithful3.pdf',width=3,height=2)
pars(mar=c(2.5,2.3,0.5,0.3))
CLT(3,dat=faithful[,'eruptions'],xlab='total duration (min)')
dev.off()

cairo(file='../../figures/CLTfaithful10.pdf',width=3,height=2)
pars()
CLT(10,dat=faithful[,'eruptions'],xlab='total duration (min)')
dev.off()

cairo(file='../../figures/galtonsbeanmachine.pdf',width=6,height=5)
pars(mar=c(2,0,0,0))
n <- 8
h <- 2
plot(c(0,n+1),c(0,n+1+h),type='n',bty='n',
     yaxt='n',xlab='',ylab='',xaxt='n')
axis(side=1,at=0.5+(0:n),labels=0:n)
for (i in 1:n){
    points((1+0.5*(i-1)):(n-0.5*(i-1)),y=h+rep(i,n-i+1))
}
mu <- (n+1)/2
minp <- dbinom(0,n,0.5)
ns <- 1/minp
k <- rep(0,n+1)
colnum <- c(4,3,5,2,6,1,7)
nc <- length(colnum)
nr <- 10
for (i in 0:n){
    p <- dbinom(i,n,0.5)
    k[i+1] <- p*ns
    lines(rep(i,2),c(0,h))
    j <- 1
    for (rn in 1:nr){
        for (cn in colnum){
            if (j>k[i+1]) break
            x <- i+cn/(nc+1)
            y <- h*rn/nr
            points(x,y,cex=0.7)
            j <- j+1
        }
        if (j>k[i+1]) break
    }
}
lines(rep(n+1,2),c(0,h))
dev.off()

draw.circle <- function(mx=.5,my=.5,radius=1,...){
    rads <- seq(from=0,to=2*pi,length.out=100)
    lines(mx+radius*sin(rads)/2,my+radius*cos(rads)/2,...)
}

get.ellipse <- function(mu=c(0.5,0.5),a=0.2,b=0.66,theta=-pi/4){
    angle <- seq(from=0,to=2*pi,length.out=100)
    rx <- mu[1] + a*cos(angle)*cos(theta) - b*sin(angle)*sin(theta)
    ry <- mu[2] + a*cos(angle)*sin(theta) + b*sin(angle)*cos(theta)
    cbind(rx,ry)
}

draw.ellipse <- function(mu=c(0.5,0.5),a=0.2,b=0.66,theta=-pi/4){
    rxry <- get.ellipse(mu=mu,a=a,b=b,theta=theta)
    lines(rxry)
}

setup_plot <- function(pop){
    plot(c(0,1),c(0,1),ann=FALSE,bty='n',type='n',asp=1)
    axis(side=1); mtext(bquote(paste('X'[.(pop)])),side=1,line=1.7,cex=0.7)
    axis(side=2); mtext(bquote(paste('Y'[.(pop)])),side=2,line=1.2,cex=0.7)    
}

cairo(file='../../figures/pop2d.pdf',width=8,height=2)
pars(mfrow=c(1,4))
setup_plot(1)
draw.circle()
setup_plot(2)
lines(c(0,1),c(0,1))
lines(c(0.5,1),c(1,1))
lines(c(1,1),c(0.5,1))
setup_plot(3)
polygon(c(0,0,1,1,0),c(0,1,1,0,0),col='black')
setup_plot(4)
draw.ellipse()
dev.off()

plot2dpop <- function(pop,n=100){
    xy <- randy(pop=pop,n=n)
    plot(xy,type='p',ann=FALSE,pch=19,cex=0.3,bty='n',xlim=c(0,1),ylim=c(0,1))
    axis(side=1); mtext(bquote(paste('X'[.(pop)])),side=1,line=1.7,cex=0.7)
    axis(side=2); mtext(bquote(paste('Y'[.(pop)])),side=2,line=1.2,cex=0.7)
}

cairo(file='../../figures/rand2d.pdf',width=8,height=2)
pars(mfrow=c(1,4))
set.seed(2)
for (i in 1:4){
    plot2dpop(i,n=100)
}
dev.off()

# returns ns sums of n samples from population pop:
plotSums <- function(n=100,ns=200,pop=1,plot=TRUE,...){
    XY <- matrix(0,nrow=ns,ncol=2) # initialise the output matrix
    for (i in 1:ns){ # compute ns sums of n values
        set.seed(i)
        xy <- randy(pop=pop,n=n) # collect n random samples
        XY[i,] <- colSums(xy) # take the sum of the X- and Y-values
    }
    if (plot){
        plot(XY,type='p',ann=FALSE,...)
        axis(side=1)
        mtext(bquote(paste('sum(X'[.(pop)]*')')),side=1,line=1.7,cex=0.7)
        axis(side=2)
        mtext(bquote(paste('sum(Y'[.(pop)]*')')),side=2,line=1.2,cex=0.7)
    }
    invisible(XY)
}

cairo(file='../../figures/rand2dsum.pdf',width=8,height=2)
pars(mfrow=c(1,4))
for (i in 1:4){
    xy <- plotSums(n=100,ns=200,pop=i,pch=19,cex=0.3,bty='n')
}
dev.off()

norm2d <- function(x,y,xym,covmat){
    out <- 0*x
    for (i in 1:length(x)){
        out[i] <- exp(-0.5*mahalanobis(c(x[i],y[i]),xym,covmat))/
            (2*pi*sqrt(det(covmat)))
    }
    out
}

cairo(file='../../figures/norm2dmarginal.pdf',width=5,height=5)
pars(mar=rep(0,4))
m <- rbind(c(1,2,2,3),c(1,4,5,3),c(1,6,7,3),c(1,8,8,3))
layout(m,widths=c(0.09,0.65,0.25,0.01),heights=c(0.01,0.25,0.65,0.09))
plot.new()
plot.new()
plot.new()
xym <- colMeans(xy)
covmat <- cov(xy)
sx <- sqrt(covmat[1,1])
sy <- sqrt(covmat[2,2])
xlim <- xym[1]+3*sx*c(-1,1)
ylim <- xym[2]+3*sy*c(-1,1)
nt <- 100
x <- seq(from=xlim[1],to=xlim[2],length.out=nt)
y <- seq(from=ylim[1],to=ylim[2],length.out=nt)
xdens <- dnorm(x,mean=xym[1],sd=sx)
plot(x,xdens,main='',xaxt='n',ann=FALSE,xlim=xlim,type='l',cex.axis=1.2)
rug(xy[,1])
mtext('Density',side=2,line=2)
plot.new()
plot(xy[,1],xy[,2],col='grey',xlab=expression('X'[4]),
     ylab=expression('Y'[4]),xlim=xlim,ylim=ylim,cex.axis=1.2)
xygrid <- expand.grid(x,y)
z <- matrix(norm2d(x=xygrid[,1],y=xygrid[,2],xym=xym,covmat=covmat),nt,nt)
contour(x,y,z,add=TRUE,labcex=0.8,xlim=xlim,ylim=ylim,
        n=100,levels=c(5,10,20,50,100,200,500)/10000,drawlabels=FALSE)
mtext(expression('X'[4]),side=1,line=2)
mtext(expression('Y'[4]),side=2,line=2)
ydens <- dnorm(y,mean=xym[2],sd=sy)
plot(ydens,y,main='',yaxt='n',ann=FALSE,ylim=ylim,type='l',cex.axis=1.2)
rug(xy[,2],side=2)
mtext('Density',side=1,line=2)
plot.new()
dev.off()

plotnormpdf <- function(mu,sigma,m,M,ylim){
    x <- seq(from=m,to=M,length.out=200)
    xlim <- c(m,M)
    plot(x,dnorm(x,mean=mu,sd=sigma),type='l',ylab='f(x)',
         xlim=xlim,ylim=ylim,bty='n')
    mtext(bquote(mu*'='*.(mu)),side=3,line=-2,at=m,adj=0,cex=0.9)
    mtext(bquote(sigma*'='*.(sigma)),side=3,line=-3,at=m,adj=0,cex=0.9)
}

cairo(file='../../figures/musigma.pdf',width=8,height=1.8)
pars(mfrow=c(1,5))
plotnormpdf(mu=0,sigma=1,m=-4,M=4,ylim=c(0,0.8))
plotnormpdf(mu=-1,sigma=1,m=-4,M=4,ylim=c(0,0.8))
plotnormpdf(mu=1,sigma=1,m=-4,M=4,ylim=c(0,0.8))
plotnormpdf(mu=0,sigma=2,m=-4,M=4,ylim=c(0,0.8))
plotnormpdf(mu=0,sigma=1/2,m=-4,M=4,ylim=c(0,0.8))
dev.off()

cairo(file='../../figures/2sigma.pdf',width=7,height=2)
pars()
layout(matrix(c(1,2,3),1,3),widths=c(0.46,0.08,0.46),heights=1)
x <- seq(from=-3,to=3,length.out=50)
plot(x,dnorm(x),type='l',bty='n',ylab='f(x)',xaxt='n',ann=FALSE)
axis(side=1,at=c(-3,-2,-1,0,1,2,3),
     labels=c(expression(mu-3*sigma),
              expression(mu-2*sigma),
              expression(mu-sigma),
              expression(mu),
              expression(mu+sigma),
              expression(mu+2*sigma),
              expression(mu+3*sigma)
              )
     )
lines(rep(-2,2),c(0,dnorm(-2)),lty=3)
lines(rep(-1,2),c(0,dnorm(-1)),lty=3)
lines(rep(0,2),c(0,dnorm(0)),lty=2)
lines(rep(1,2),c(0,dnorm(1)),lty=3)
lines(rep(2,2),c(0,dnorm(2)),lty=3)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(1,0))
p <- par(mar=c(2.5,0,0.5,0))
plot(c(0,1),c(0,1),type='n',ann=FALSE,xaxt='n',yaxt='n',bty='n',xpd=TRUE)
lines(rep(0.5,2),pnorm(c(-2,2)))
text(0.35,0.5,labels='95%',srt=90,pos=NULL)
lines(rep(1,2),pnorm(c(-1,1)))
text(0.85,0.5,labels='68%',srt=90,pos=NULL)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
par(mar=c(2.5,2.3,0.5,0.5))
plot(x,pnorm(x),type='l',ylim=c(0,1),bty='n',
     ann=FALSE,xaxt='n',yaxt='n',xpd=TRUE)
axis(side=1,at=c(-3,-2,-1,0,1,2,3),
     labels=c(expression(mu-3*sigma),
              expression(mu-2*sigma),
              expression(mu-sigma),
              expression(mu),
              expression(mu+sigma),
              expression(mu+2*sigma),
              expression(mu+3*sigma))
     )
lines(c(rep(-2,2),-3),c(0,rep(pnorm(-2),2)),lty=3)
lines(c(rep(-1,2),-3),c(0,rep(pnorm(-1),2)),lty=3)
lines(c(rep(0,2),-3),c(0,rep(pnorm(0),2)),lty=2)
lines(c(rep(1,2),-3),c(0,rep(pnorm(1),2)),lty=3)
lines(c(rep(2,2),-3),c(0,rep(pnorm(2),2)),lty=3)
axis(side=2,at=pnorm(c(-2,-1,0,1,2)),
     labels=c(signif(pnorm(-2),3),
              signif(pnorm(-1),3),
              signif(pnorm(0),3),
              signif(pnorm(1),3),
              signif(pnorm(2),3)),
     las=2
     )
lines(rep(-4,2),c(0,1))
par(p)
dev.off()

plot2dnorm <- function(xym,covmat,...){
    nt <- 100
    m <- -5
    M <- 5
    x <- seq(from=m,to=M,length=nt)
    y <- seq(from=m,to=M,length=nt)
    xygrid <- expand.grid(x,y)
    z <- norm2d(x=xygrid[,1],y=xygrid[,2],xym=xym,covmat=covmat)
    contour(x=x,y=y,z=matrix(z,nt,nt),labcex=0.8,xlab='X',
            asp=1,bty='n',drawlabels=FALSE,nlevels=6,...)
    lines(c(-1,1),rep(m,2))
    sxlab <- bquote(sigma[x]*'='*.(sqrt(covmat[1,1])))
    text(xym[1],m,labels=sxlab,pos=3,offset=0.1)
    lines(rep(m,2),c(-1,1))
    sylab <- bquote(sigma[y]*'='*.(sqrt(covmat[2,2])))
    text(m,xym[2],labels=sylab,pos=4,offset=0.1)
    sxylab <- bquote(sigma[x*','*y]*'='*.(covmat[1,2]))
    legend('topleft',legend=sxylab,bty='n')
}

cairo(file='../../figures/cov.pdf',width=7,height=1.75)
pars(mfrow=c(1,4))
plot2dnorm(xym=c(0,0),covmat=matrix(c(1,0,0,1),2,2),ylab='Y')
plot2dnorm(xym=c(0,0),covmat=matrix(c(1/4,0,0,4),2,2),yaxt='n')
plot2dnorm(xym=c(0,0),covmat=matrix(c(1,.9,.9,1),2,2),yaxt='n')
plot2dnorm(xym=c(0,0),covmat=matrix(c(1,-.9,-.9,1),2,2),yaxt='n')
dev.off()

darts <- function(mu,covmat,n=10){
    col <- 'gray50'
    plot(c(-1,1),c(-1,1),type='n',bty='n',ann=FALSE,asp=1,xaxt='n',yaxt='n')
    draw.circle(mx=0,my=0,radius=0.2,col=col)
    draw.circle(mx=0,my=0,radius=1,col=col)
    draw.circle(mx=0,my=0,radius=2,col=col)
    xy <- mvrnorm(n,mu=mu,Sigma=covmat)
    points(xy[,1],xy[,2],pch=21,cex=0.8,lwd=1.5,col='black',bg='white')
}

cairo(file='../../figures/accuracyvsprecision.pdf',width=7,height=2)
pars(mar=rep(0,4))
m <- rbind(c(1,1,1,1),c(2,3,4,5),c(6,6,6,6))
layout(m,widths=1,heights=c(0.2,0.6,0.2))
plot.new()
set.seed(7)
cex <- 0.9
darts(mu=c(0,0),covmat=.01*matrix(c(1,0,0,1),2,2),n=20)
mtext('high precision',side=3,line=1.5,cex=cex)
mtext('high accuracy',side=3,line=0,cex=cex)
darts(mu=c(0,0),covmat=.05*matrix(c(1,0,0,1),2,2),n=20)
mtext('low precision',side=3,line=1.5,cex=cex)
mtext('high accuracy',side=3,line=0,cex=cex)
darts(mu=c(.5,.5),covmat=.05*matrix(c(1,0,0,1),2,2),n=20)
mtext('low precision',side=3,line=1.5,cex=cex)
mtext('low accuracy',side=3,line=0,cex=cex)
darts(mu=c(.5,.5),covmat=.01*matrix(c(1,0,0,1),2,2),n=20)
mtext('high precision',side=3,line=1.5,cex=cex)
mtext('low accuracy',side=3,line=0,cex=cex)
plot(c(0,1),c(0,1),type='n',bty='n',xaxt='n',yaxt='n',ann=FALSE)
arrows(0.1,.5,0.9,.5)
text(0,0.5,labels='good',cex=1.4)
text(1,0.5,labels='bad',cex=1.4)
dev.off()

cairo(file='../../figures/errorprop2d.pdf',width=3,height=3)
pars(mar=rep(0,4))
m <- rbind(rep(7,4),c(1,2,5,8),c(1,3,4,8),c(c(1,6,6,8)))
layout(m,widths=c(0.11,0.7,0.2,0.02),heights=c(0.02,0.2,0.7,0.11))
plot.new()
nt <- 200
m <- -10
M <- 10
x <- seq(from=m,to=5,length.out=nt)
y <- dnorm(x,-5)
xlim <- c(m,5)
ylim <- c(-5,M)
plot(x,y,type='l',bty='n',xaxt='n',yaxt='n',ann=FALSE)
lines(x,0*x)
axis(side=2,at=c(0,0.1,0.2,0.3,0.4),labels=c('0',' ','0.2',' ','0.4'))
mtext('f(x)',side=2,line=1.2,cex=0.7)
f1 <- function(x){5-x/2}
x1 <- seq(from=xlim[1],to=xlim[2],length.out=nt)
y1 <- f1(x1)
plot(x1,y1,type='l',bty='n',xaxt='n',yaxt='n',ann=FALSE,xlim=xlim,ylim=ylim)
f2 <- function(x){-10-2*x}
x2 <- seq(from=m,to=-2.5,length.out=nt)
y2 <- f2(x2)
lines(x2,y2,col='gray50')
lines(rep(-5,2),c(9.5,f1(-5)),lty=2,col='black')
lines(rep(-3,2),c(9.5,f1(-3)),lty=2,col='black')
lines(rep(-5,2),c(f1(-5),f2(-5)),lty=2,col='gray50')
lines(rep(-3,2),c(f1(-3),f2(-3)),lty=2,col='gray50')
lines(c(-5,4),rep(f1(-5),2),lty=2)
lines(c(-3,4),rep(f1(-3),2),lty=2)
lines(c(-5,4),rep(f2(-5),2),lty=2,col='gray50')
lines(c(-3,4),rep(f2(-3),2),lty=2,col='gray50')
text(-5,10,expression(bar('x')))
text(-3,10,expression('x'[i]))
text(5,7.5,expression(bar('z')[1]))
text(5,6.5,expression('z'[1*i]))
text(5,0,expression(bar('z')[2]),col='gray50')
text(5,-4,expression('z'[2*i]),col='gray50')
text(-7.5,f1(-7)+.2,expression('g'[1]*'(x)'),pos=1)
text(-7,f2(-7),expression('g'[2]*'(x)'),col='gray50',pos=2)
axis(side=1)
mtext('x',side=1,line=1.3,cex=0.7)
axis(side=2)
mtext('z',side=2,line=1.2,cex=0.7)
yy <- seq(from=-5,to=10,length.out=nt)
d1 <- dnorm(yy,f1(-5),1/2)
d2 <- dnorm(yy,f2(-5),2)
plot(d1,yy,type='l',yaxt='n',bty='n',xaxt='n',ann=FALSE,ylim=ylim)
lines(d2,yy,col='gray50')
lines(0*yy,yy)
axis(side=1,at=c(0,0.2,0.4,0.6,0.8),labels=c('0',' ','0.4',' ','0.8'))
mtext('f(z)',side=1,line=1.3,cex=0.7)
plot.new()
plot.new()
dev.off()

cairo(file='../../figures/errorprop3d.pdf',width=5.5,height=5.5)
pars()
nt <- 30
xm <- -0.2
ym <- -0.8
m <- -1
M <- 1
x <- seq(from=m,to=M,length.out=nt)
y <- seq(from=m,t=M,length.out=nt)
xy <- expand.grid(x,y)
pol <- function(xy){
    xy[,1]*xy[,2]^2
}
z <- pol(xy)
p <- persp(x,y,matrix(z,nt,nt),phi=30,theta=-30,
           zlab='z',border='gray50')
ell <- get.ellipse(mu=c(xm,ym),a=0.07,b=0.15,theta=-pi/3)
#
xyell1 <- trans3d(x=ell[,1],y=ell[,2],z=rep(-1,nrow(ell)),pmat=p)
lines(xyell1$x,xyell1$y)
xyp1 <- trans3d(x=xm,y=ym,z=-1,pmat=p)
points(xyp1$x,xyp1$y,pch=19,cex=0.5)
#
zell <- pol(ell)
xyell2 <- trans3d(x=ell[,1],y=ell[,2],z=zell,pmat=p)
lines(xyell2$x,xyell2$y)
z <- pol(rbind(c(xm,ym)))
xyp2 <- trans3d(x=xm,y=ym,z=z,pmat=p)
points(xyp2$x,xyp2$y,pch=19,cex=0.5)
#
xyell3 <- trans3d(x=M,y=ell[,2],z=zell,pmat=p)
lines(xyell3$x,xyell3$y)
xyp3 <- trans3d(x=M,y=ym,z=z,pmat=p)
points(xyp3$x,xyp3$y,pch=19,cex=0.5)
#
xm1 <- which.min(xyell1$x)
xM1 <- which.max(xyell1$x)
ym1 <- which.min(xyell1$y)
yM1 <- which.max(xyell1$y)
xm2 <- which.min(xyell2$x)
xM2 <- which.max(xyell2$x)
ym2 <- which.min(xyell2$y)
yM2 <- which.max(xyell2$y)
ym3 <- which.min(xyell3$y)
yM3 <- which.max(xyell3$y)
lines(x=c(xyell1$x[xm1],xyell2$x[xm2]),
      y=c(xyell1$y[xm1],xyell2$y[xm2]),lty=3)
lines(x=c(xyell1$x[xM1],xyell2$x[xM2]),
      y=c(xyell1$y[xM1],xyell2$y[xM2]),lty=3)
lines(x=c(xyell2$x[ym2],xyell3$x[ym3]),
      y=c(xyell2$y[ym2],xyell3$y[ym3]),lty=3)
lines(x=c(xyell2$x[yM2],xyell3$x[yM3]),
      y=c(xyell2$y[yM2],xyell3$y[yM3]),lty=3)
#
drx <- 10
mxax <- trans3d(x=ell[xm1+drx,1],y=m,z=-1,pmat=p)
Mxax <- trans3d(x=ell[xM1+drx,1],y=m,z=-1,pmat=p)
lines(x=c(mxax$x,xyell1$x[xm1+drx]),y=c(mxax$y,xyell1$y[xm1+drx]),lty=2)
lines(x=c(Mxax$x,xyell1$x[xM1+drx]),y=c(Mxax$y,xyell1$y[xM1+drx]),lty=2)
dry <- 10
myax <- trans3d(x=M,y=ell[ym3+dry,2],z=-1,pmat=p)
Myax <- trans3d(x=M,y=ell[yM3+dry,2],z=-1,pmat=p)
lines(x=c(myax$x,xyell3$x[ym3+dry]),y=c(myax$y,xyell3$y[ym3+dry]),lty=2)
lines(x=c(Myax$x,xyell3$x[yM3+dry]),y=c(Myax$y,xyell3$y[yM3+dry]),lty=2)
drxy <- 10
lines(x=c(myax$x,xyell1$x[ym1+drxy]),y=c(myax$y,xyell1$y[ym1+drxy]),lty=2)
lines(x=c(Myax$x,xyell1$x[yM1+drxy]),y=c(Myax$y,xyell1$y[yM1+drxy]),lty=2)
drz <- -5
zm <- which.min(zell)
zM <- which.max(zell)
mzax <- trans3d(x=M,y=m,z=zell[zm+drz],pmat=p)
Mzax <- trans3d(x=M,y=m,z=zell[zM+drz],pmat=p)
lines(x=c(mzax$x,xyell3$x[ym3+drz]),y=c(mzax$y,xyell3$y[ym3+drz]),lty=2)
lines(x=c(Mzax$x,xyell3$x[yM3+drz]),y=c(Mzax$y,xyell3$y[yM3+drz]),lty=2)
#
txy1 <- trans3d(x=xm,y=m,z=-1,pmat=p)
txy2 <- trans3d(x=M,y=ym,z=-1,pmat=p)
txy3 <- trans3d(x=M,y=m,z=z,pmat=p)
text(txy1$x-0.002,txy1$y+0.013,labels=expression(s[x]),pos=1,srt=30)
text(txy2$x+0.005,txy2$y+0.01,labels=expression(s[y]),pos=3,srt=-20)
text(txy3$x,txy3$y-0.01,labels=expression(s[z]),pos=4,offset=0.1,srt=-20)
dev.off()

if (FALSE){ # Q-Q plot of earthquake data
    nvals <- length(unique(nquakes))
    qqplot(qpois(ppoints(nvals),lambda=mean(nquakes)),
           nquakes)
    qqline(nquakes,
           distribution = function(probs) { qpois(probs, lambda=lambda) },
           col = "red",
           lwd = 0.5)
}

qqfaithful <- function(n,fname){
    cairo(file=fname,width=3,height=2.5)
    pars()
    if (n==1) clt <- faithful[,'eruptions']
    else clt <- CLT(n,dat=faithful[,'eruptions'],plot=FALSE)
    qqnorm(clt,main='',cex=0.8)
    qqline(clt)
    dev.off()
}
qqfaithful(n=1,'../../figures/qqfaithful1.pdf')
qqfaithful(n=2,'../../figures/qqfaithful2.pdf')
qqfaithful(n=3,'../../figures/qqfaithful3.pdf')
qqfaithful(n=10,'../../figures/qqfaithful10.pdf')

cairo(file='../../figures/qqfaithful12.pdf',width=3,height=2.5)
pars(mar=c(2.5,2.5,0.5,0.25))
xy1 <- plotSums(n=100,ns=200,pop=1,plot=FALSE)
xy2 <- plotSums(n=100,ns=200,pop=2,plot=FALSE)
qqplot(xy1[,1],xy2[,1],
       xlab=expression('x'[1]*'=sum(X'[1]*')'),
       ylab=expression('x'[2]*'=sum(X'[2]*')'))
dev.off()

qtiles <- c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99)
#signif(qt(qtiles,df=4),2)

cairo(file='../../figures/1samplettest.pdf',width=6,height=3)
set.seed(1)
n <- 5
df <- n-1
gold1 <- signif(sort(rnorm(n,19.15,0.1)),df)
tt <- t.test(gold1,mu=19.3)
pars(mfrow=c(1,2))
x <- seq(from=qt(0.001,df=tt$parameter),
         to=qt(0.999,df=tt$parameter),length.out=100)
y <- dt(x,df=tt$parameter)
plot(x,y,type='l',bty='n',xlab='t',ylab='f(t)')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1))
lines(rep(tt$statistic,2),range(y),lty=2)
ltail <- (x<qt(0.025,df=df))
polygon(c(x[ltail][1],x[ltail],tail(x[ltail],1)),c(0,y[ltail],0),col='black')
Y <- pt(x,df=df)
plot(x,Y,type='l',bty='n',xlab='t',ylab='F(t)')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(1,1))
lines(rep(tt$statistic,2),range(Y),lty=2)
lines(range(x),rep(pt(tt$statistic,df=df),2),lty=2)
lines(range(x),rep(0.05,2),lty=3)
dev.off()

gold2 <- signif(sort(rnorm(n-1,19.25,0.1)),df)

cairo(file='../../figures/2samplettest.pdf',width=6,height=3)
set.seed(1)
df <- 7
gold1 <- signif(sort(rnorm(n,19.15,0.1)),df)
#signif(qt(qtiles,df=df),2)
tt <- t.test(gold1,gold2)
pars(mfrow=c(1,2))
x <- seq(from=qt(0.001,df=tt$parameter),
         to=qt(0.999,df=tt$parameter),length.out=100)
y <- dt(x,df=df)
plot(x,y,type='l',bty='n',xlab='t',ylab='f(t)')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1))
lines(rep(tt$statistic,2),range(y),lty=2)
ltail <- (x<qt(0.025,df=df))
polygon(c(x[ltail][1],x[ltail],tail(x[ltail],1)),c(0,y[ltail],0),col='black')
utail <- (x>qt(0.975,df=df))
polygon(c(x[utail][1],x[utail],tail(x[utail],1)),c(0,y[utail],0),col='black')
Y <- pt(x,df=df)
plot(x,Y,type='l',bty='n',xlab='t',ylab='F(t)')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(1,1))
lines(rep(tt$statistic,2),range(Y),lty=2)
lines(range(x),rep(pt(tt$statistic,df=df),2),lty=2)
lines(range(x),rep(0.025,2),lty=3)
lines(range(x),rep(0.975,2),lty=3)
dev.off()

cairo(file='../../figures/normconf.pdf',width=6,height=3)
pars(mfrow=c(1,2))
x <- seq(from=-5,to=5,length.out=100)
yl <- dnorm(x,mean=-2)
yu <- dnorm(x,mean=2)
plot(x,yl,type='l',bty='n',col='gray50',xaxt='n',xlab='',ylab='f(x)')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(1,1))
lines(x,yu,type='l')
lines(rep(0,2),range(yl),lty=2)
lines(rep(-2,2),range(yl),lty=3,col='gray50')
lines(rep(2,2),range(yl),lty=3)
axis(side=1,at=c(-4,-2,0,2,4),
     labels=c('',
              expression(bar('x')*'-2s['*bar('x')*']'),
              expression(bar('x')),
              expression(bar('x')*'+2s['*bar('x')*']'),
              '')
     )
Yl <- pnorm(x,mean=-2)
Yu <- pnorm(x,mean=2)
plot(x,Yl,type='l',bty='n',col='gray50',xaxt='n',xlab='',ylab='F(x)')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(1,1))
lines(x,Yu,type='l')
lines(rep(0,2),range(Yl),lty=2)
lines(rep(-2,2),range(Yl),lty=3,col='gray50')
lines(rep(2,2),range(Yu),lty=3)
lines(range(x),rep(0.025,2),lty=3)
lines(range(x),rep(0.975,2),lty=3)
axis(side=1,at=c(-4,-2,0,2,4),
     labels=c('',
              expression(bar('x')*'-2s['*bar('x')*']'),
              expression(bar('x')),
              expression(bar('x')*'+2s['*bar('x')*']'),
              '')
     )
dev.off()

cairo(file='../../figures/tdof.pdf',width=6,height=3)
pars(mfrow=c(1,2))
x <- seq(from=qt(0.1,df=1),
         to=qt(0.9,df=1),length.out=100)
plot(range(x),c(0,0.4),type='n',bty='n',xlab='t',ylab='f(t)')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(1,0))
legend('bottom',legend=c('df=1','df=3','df=30'), lty=c(3,2,1),
       cex=1.0,inset=c(0.05,0.05),bty='n')
lines(x,dt(x,df=1),lty=3)
lines(x,dt(x,df=3),lty=2)
lines(x,dt(x,df=30),lty=1)
plot(range(x),c(0,1),type='n',bty='n',xlab='t',ylab='F(t)')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(1,0))
lines(x,pt(x,df=1),lty=3)
lines(x,pt(x,df=3),lty=2)
lines(x,pt(x,df=30),lty=1)
legend('bottomright',legend=c('df=1','df=3','df=30'),
       lty=c(3,2,1),cex=1.0,inset=c(0.05,0.05),bty='n')
dev.off()

# signif(-qt(0.025,df=c(1:10,30,100,1000)),4)

quaketab <- table(nquakes)
pp <- dpois(0:25,lambda=mean(nquakes))
obs <- c(sum(quaketab[1:2]),quaketab[3:9],sum(quaketab[10:12]))
pred <- c(sum(pp[1:3]),pp[4:10],sum(pp[11:20]))
p <- pred/sum(pred)
o <- obs/sum(obs)
X <- chisq.test(obs,p=p)
w <- sqrt(sum((o-p)^2/p))

# signif(qchisq(qtiles,df=5),2)

cairo(file='../../figures/chi2.pdf',width=6,height=3)
pars(mfrow=c(1,2))
x <- seq(from=0,to=25,length.out=100)
y <- dchisq(x,df=X$parameter-1)
xmax <- qchisq(0.95,df=X$parameter-1)
plot(x,y,type='l',bty='n',xlab=expression(chi^2),ylab='')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,-.6))
mtext(expression('f['*chi^2*']'),side=2,line=1.1)
utail <- x>xmax
polygon(c(x[utail][1],x[utail],tail(x[utail],1)),
        c(0,y[utail],0),col='black')
lines(rep(X$statistic,2),range(y),lty=2)
y <- pchisq(x,df=X$parameter-1)
plot(x,y,type='l',bty='n',xlab=expression(chi^2),ylab='')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,-.6))
mtext(expression('F['*chi^2*']'),side=2,line=1.1)
lines(rep(X$statistic,2),range(y),lty=2)
lines(range(x),rep(0.95,2),lty=3)
pval <- pchisq(X$statistic,df=X$parameter-1)
lines(range(x),rep(pval,2),lty=2)
dev.off()

zirctab <- table(nzirc)
obs <- c(zirctab[1:5],sum(zirctab[6:8]))
pp <- dpois(0:20,lambda=mean(nzirc))
pred <- c(sum(pp[1:2]),pp[3:6],sum(pp[7:20]))
p <- pred/sum(pred)
X <- chisq.test(obs,p=p)
xmax <- qchisq(0.95,df=X$parameter)

cairo(file='../../figures/cherrypicking.pdf',width=10,height=5)
pars(mar=rep(0,4))
m <- rbind(c(0.04,0.23,0.60,0.95),
           c(0.27,0.48,0.60,0.95),
           c(0.52,0.73,0.60,0.95),
           c(0.77,0.96,0.60,0.95),
           c(0.20,0.80,0.1,0.52),
           c(0.00,1.00,0.00,1.00))
split.screen(m)
screen(1); par(mar=rep(0,4))
d1 <- dpois(1:6,lambda=mean(nzirc))
b <- barplot(48*d1,xlab='n',ylab='n',xaxt='n',yaxt='n',main='',col='grey90')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
axis(side=1,at=b,labels=c('<2',2,3,4,5,'>5'))
axis(side=2)
#mtext(paste0('predicted distribution'),side=3,line=0.1)
screen(2,new=FALSE); par(mar=rep(0,4))
b <- barplot(48*d1,xlab='n',ylab='n',xaxt='n',yaxt='n',main='',col='white')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
axis(side=1,at=b,labels=c('<2',2,3,4,5,'>5'))
axis(side=2)
mtext(expression(chi^2*'=0, p=1.0'),side=3,line=-0.1)
screen(3,new=FALSE); par(mar=rep(0,4))
barplot(obs,xlab='n',ylab='n',xaxt='n',yaxt='n',main='',col='white')
legend('topleft',legend='c)',bty='n',cex=1.2,adj=c(2,0))
axis(side=1,at=b,labels=c('<2',2,3,4,5,'>5'))
axis(side=2)
X1 <- chisq.test(obs,p=p)$statistic
mtext(bquote(chi^2*'='*.(signif(X1,2))*', p='*.(signif(1-pchisq(X1,df),2))),
      side=3,line=-0.1)
screen(4,new=FALSE); par(mar=rep(0,4))
barplot(c(3,2,3,12,20,8),xlab='n',ylab='n',xaxt='n',yaxt='n',main='',col='white')
legend('topleft',legend='d)',bty='n',cex=1.2,adj=c(2,0))
axis(side=2)
axis(side=1,at=b,labels=c('<2',2,3,4,5,'>5'))
X2 <- chisq.test(c(3,2,3,12,20,8),p=p)$statistic
mtext(bquote(chi^2*'='*.(signif(X2,2))*', p='*.(signif(1-pchisq(X2,df),2))),
      side=3,line=-0.1)
screen(5,new=FALSE); par(mar=rep(0,4))
x <- seq(from=0,to=50,length.out=100)
y <- dchisq(x,df=5)
plot(x,y,type='l',bty='n',xlab='',ylab='')
legend('topleft',legend='e)',bty='n',cex=1.2,adj=c(2,0))
mtext(expression(chi^2),side=1,line=1.5)
mtext(expression('p('*chi^2*')'),side=2,line=2)
lines(rep(X1,2),c(0,dchisq(X1,df=5)),lt=2)
utail <- x>xmax
polygon(c(x[utail][1],x[utail],tail(x[utail],1)),
        c(0,y[utail],0),col='black')
screen(6,new=FALSE); par(mar=rep(0,4))
plot(c(0,1),c(0,1),type='n',bty='n',ann=FALSE,xaxt='n',yaxt='n')
lines(x=c(0.20,0.36),y=c(0.083,0.53),lty=2)
lines(x=c(0.290,0.64),y=c(0.268,0.53),lty=2)
lines(x=c(0.745,0.90),y=c(0.083,0.53),lty=2)
close.screen(6,all.screens=TRUE)
dev.off()

A <- c(10,5,6,20)
B <- c(25,12,10,35)
clasts <- rbind(A,B)
p <- prop.table(clasts)
rs <- rowSums(clasts)
cs <- colSums(clasts)
expected <- signif(matrix(rs,nrow=2)%*%matrix(cs,ncol=4)/sum(clasts),3)
X <- chisq.test(clasts,expected/sum(expected))

# signif(qchisq(qtiles,df=3),2)

cairo(file='../../figures/chi22.pdf',width=6,height=3)
pars(mfrow=c(1,2),mar=c(2.5,2.3,0.6,0.2))
x <- seq(from=0,to=15,length.out=100)
y <- dchisq(x,df=X$parameter)
xmax <- qchisq(0.95,df=X$parameter)
plot(x,y,type='l',bty='n',xlab=expression(chi^2),ylab='')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,-1.2),xpd=NA)
mtext(expression('f['*chi^2*']'),side=2,line=1.1)
utail <- x>xmax
polygon(c(x[utail][1],x[utail],tail(x[utail],1)),
        c(0,y[utail],0),col='black')
lines(rep(X$statistic,2),c(0,0.24),lty=2)
y <- pchisq(x,df=X$parameter)
plot(x,y,type='l',bty='n',xlab=expression(chi^2),ylab='')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,-1.2),xpd=NA)
mtext(expression('F['*chi^2*']'),side=2,line=1.1)
lines(rep(X$statistic,2),range(y),lty=2)
lines(range(x),rep(0.95,2),lty=3)
pval <- pchisq(X$statistic,df=3)
lines(range(x),rep(pval,2),lty=2)
dev.off()

# signif(qchisq(qtiles,df=6),3)

p <- c(0.2845528,0.1382114,0.1300813,0.4471545)
set.seed(1)
A <- round(100000*p*(1 + (0.5+runif(4))/20))
B <- round(100000*p*(1 + (0.5+runif(4))/20))
AB <- rbind(A,B)           
rs <- rowSums(AB)
cs <- colSums(AB)
expected <- signif(matrix(rs,nrow=nrow(AB))%*%matrix(cs,ncol=ncol(AB))/sum(AB),5)
X <- chisq.test(AB,expected/sum(expected))

# signif(9+qwilcox(qtiles,n=5,m=4),3)

cairo(file='../../figures/wilcox.pdf',width=6,height=3)
pars(mfrow=c(1,2))
x <- 9:36
y <- dwilcox(x-10,m=5,n=4)
llim <- qwilcox(0.025,4,5)+1
ulim <- qwilcox(0.975,4,5)+3
col <- rep('white',length(x))
col[1:llim] <- 'black'
col[ulim:length(x)] <- 'black'
b <- barplot(y,xaxt='n',col=col,space=0,ylab='P(W)')
lines(rep(26,2)-8.5,range(y),lty=2)
ticks <- c(9,15,20,25,30,36)
axis(side=1,
     at=b[1,1]+(b[nrow(b),1]-b[1,1])*(ticks-9)/(36-9),
     labels=ticks)
mtext('W',side=1,line=1.5)
xx <- x
yy <- pwilcox(xx-10,m=5,n=4)
plot(xx,yy,type='s',bty='n',xaxt='n',xlab='x',ylab='F(x)')
lines(rep(26,2),range(yy),lty=2)
lines(range(xx),rep(0.025,2),lty=3)
lines(range(xx),rep(0.975,2),lty=3)
lines(range(xx),rep(pwilcox(25-9,n=5,m=4),2),lty=2)
axis(side=1,at=ticks)
dev.off()

cairo(file='../../figures/KS.pdf',width=3.5,height=3.5)
pars()
na.rm <- function(dat,var){
    dat[!is.na(dat[,var]),var]
}
DZ <- read.csv('DZages.csv',header=TRUE)
samp1 <- na.rm(DZ,'X5')
samp2 <- na.rm(DZ,'Y')
cdf1 <- ecdf(samp1)
cdf2 <- ecdf(samp2)
plot(cdf2,verticals=TRUE,pch=NA,main='',xlab='age (Ma)',ylab='F(age)')
plot(cdf1,verticals=TRUE,pch=NA,add=TRUE,col='gray50')
x1 <- sort(samp1)
y1 <- seq(from=0,to=1,length.out=length(x1))
D <- 0
for (i in 1:length(x1)){
    d <- abs(cdf2(x1[i])-y1[i])
    if (d>D) {
        D <- d
        xi <- x1[i]
    }
}
x2 <- sort(samp2)
y2 <- seq(from=0,to=1,length.out=length(x2))
for (i in 1:length(x2)){
    d <- abs(cdf1(x2[i])-y2[i])
    if (d>D){
        D <- d
        xi <- x2[i]
    }
}
ks <- ks.test(samp1,samp2)
lines(rep(xi,2),c(cdf1(xi),cdf2(xi)),lwd=1.5)
text(xi,cdf1(xi)*.3+cdf2(xi)*.7,pos=2,srt=90,
     labels=paste0('D=',signif(ks$statistic,2)))
text(1000,cdf1(1000),pos=1,labels='dune',offset=1,col='gray50')
text(1000,cdf2(1000),pos=3,labels='river')
dev.off()

cairo(file='../../figures/KSdens.pdf',width=6,height=3)
pars(mfrow=c(1,2))
n.x <- length(samp1)
n.y <- length(samp2)
y1 <- seq(from=0,to=1,length.out=n.x)
y2 <- seq(from=0,to=1,length.out=n.y)
y12 <- expand.grid(y1,y1)
dy <- sort(abs(y12[,1]-y12[,2]))
uniek <- which(diff(dy)>1e-10)
D <- dy[uniek]
nD <- length(D)
d <- rep(0,2*nD)
P <- 0*D
p <- rep(0,2*nD)
for (i in 1:nD){
    P[i] <- .Call(stats:::C_pSmirnov2x, D[i],n.x=n.x,n.y=n.y)
    p[2*i+c(-1,0)] <- P[i]
    d[2*i+c(-1,0)] <- D[i]
}
dens <- diff(P)
rej <- which(P>0.95)
col <- rep('white',length(dens))
col[rej] <- 'black'
ticks <- c(0,0.1,0.2,0.3,0.4)
toplot <- which(D<0.4)
b <- barplot(dens[toplot],space=0,col=col)
at <- b[1,1]+ticks*(b[nrow(b),1]-b[1,1])/diff(range(ticks))
axis(side=1,at=at,labels=ticks)
mtext('D',side=1,line=1.5)
mtext('f(D)',side=2,line=1.5)
bd <- b[1,1]+ks$statistic*(b[nrow(b),1]-b[1,1])/diff(range(ticks))
lines(rep(bd,2),range(dens[toplot],2),lty=2)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1.2),xpd=NA)
plot(d[-1],p[-length(p)],type='l',xlim=range(ticks),bty='n',xlab='x',ylab='F(x)')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,1.2),xpd=NA)
lines(rep(ks$statistic,2),c(0,1),lty=2)
lines(range(d[-1]),rep(0.95,2),lty=3)
lines(range(d[-1]),rep(1,2),lty=2)
lookup <- function(q,D,P){
    rej <- which(P>q)
    D[min(rej)]
}
q95 <- lookup(0.95,D,P)
lines(rep(q95,2),c(0,1),lty=3)
dev.off()

theD <- signif(sapply(qtiles,lookup,D=D,P=P),3)

cairo(file='../../figures/KSnorm.pdf',width=3.5,height=3.5)
pars()
xlim <- c(-1500,4000)
plot(ecdf(samp1),pch=NA,verticals=TRUE,xlim=xlim,main='')
x <- seq(from=xlim[1],to=xlim[2],length.out=100)
y <- pnorm(x,mean(samp1),sd(samp1))
py <- pnorm(samp1,mean(samp1),sd(samp1))
D <- max(abs(py-rank(samp1)/length(samp1)))
Fsamp1 <- rank(samp1)/length(samp1)
i <- which.max(abs(py-Fsamp1))
lines(x,y,col='gray50')
lines(rep(samp1[i],2),c(py[i],Fsamp1[i]),lwd=1.5)
text(samp1[i],.5*py[i]+.5*Fsamp1[i],labels='D',srt=0,pos=4,offset=0.1)
dev.off()

RbSrGenerator <- function(n){
    lambda <- 1.42e-5
    tt <- 1000
    set.seed(5)
    RbSr <- 1+runif(n)*9
    rho <- 0.5*(1+rnorm(n,0,1)/10)
    e <- matrix(0,n,n)
    E <- matrix(0,2,2)
    SrSr0 <- 0.7
    SrSr <- SrSr0 + RbSr*(exp(lambda*tt)-1)
    errRbSr <- 0*RbSr
    errSrSr <- 0*SrSr
    for (i in 1:n){
        E[1,1] <- (1+rnorm(1,0,1)/10)*1e-2
        E[2,2] <- (1+rnorm(1,0,1)/10)*5e-5
        E[1,2] <- rho[i]*sqrt(E[1,1]*E[2,2])
        E[2,1] <- E[1,2]
        e <- mvrnorm(1,rep(0,2),E)
        RbSr[i] <- RbSr[i] + e[1]
        SrSr[i] <- SrSr[i] + e[2]
        errRbSr[i] <- sqrt(E[1,1])
        errSrSr[i] <- sqrt(E[2,2])
    }
    signif(rbind(RbSr,SrSr,errRbSr,errSrSr,rho),3)    
}

RbSr <- RbSrGenerator(n=8)

a <- 0.6
b <- 0.03
yhat <- a+b*RbSr['RbSr',]
e <- yhat-RbSr['SrSr',]
tab1 <- signif(rbind(RbSr['RbSr',],RbSr['SrSr',],yhat,e),3)
tab2 <- signif(rbind(tab1[1:2,],tab1[1,]^2,tab1[1,]*tab1[2,]),3)
tab2 <- cbind(tab2,rowSums(tab2))

cairo(file='../../figures/r.pdf',width=7.2,height=1.9)
pars(mfrow=c(1,4))
set.seed(2)
n <- 100
r <- 0; xy <- mvrnorm(n,rep(0,2),rbind(c(1,r),c(r,1)))
plot(xy,pch=19,asp=1,xlim=c(-3,3),ylim=c(-3,3.5),bty='n',axes=FALSE,xlab='x',ylab='y',cex=.7)
legend(1,-1,legend=c(paste0('r=',signif(r,2)),expression('r'^2*'=0')),bty='n')
mtext(side=3,text='a)',adj=0,line=-2,cex=0.8)
axis(side=1,at=c(-2,-1,0,1,2))
axis(side=2,at=c(-2,-1,0,1,2))
r <- sqrt(0.5); xy <- mvrnorm(n,rep(0,2),rbind(c(1,r),c(r,1)))
plot(xy,pch=19,asp=1,xlim=c(-3,3),ylim=c(-3,3.5),bty='n',axes=FALSE,xlab='x',ylab='',cex=.7)
legend(-0.5,-1,legend=c(paste0('r=',signif(r,2)),expression('r'^2*'=0.5')),bty='n')
mtext(side=3,text='b)',adj=0,line=-2,cex=0.8)
axis(side=1,at=c(-2,-1,0,1,2))
r <- sqrt(0.9); xy <- mvrnorm(n,rep(0,2),rbind(c(1,r),c(r,1)))
plot(xy,pch=19,asp=1,xlim=c(-3,3),ylim=c(-3,3.5),bty='n',axes=FALSE,xlab='x',ylab='',cex=.7)
mtext(side=3,text='c)',adj=0,line=-2,cex=0.8)
legend(-0.5,-1,legend=c(paste0('r=',signif(r,2)),expression('r'^2*'=0.9')),bty='n')
axis(side=1,at=c(-2,-1,0,1,2))
r <- -sqrt(0.9); xy <- mvrnorm(n,rep(0,2),rbind(c(1,r),c(r,1)))
plot(xy,pch=19,asp=1,xlim=c(-3,3),ylim=c(-3,3.5),bty='n',axes=FALSE,xlab='x',ylab='',cex=.7)
mtext(side=3,text='d)',adj=0,line=-2,cex=0.8)
legend(x=-3,y=-1,legend=c(paste0('r=',signif(r,2)),expression('r'^2*'=0.9')),bty='n')
axis(side=1,at=c(-2,-1,0,1,2))
dev.off()

r <- cor(t(RbSr[c('RbSr','SrSr'),]))[1,2]
n <- ncol(RbSr)
tt <- r*sqrt(n-2)/sqrt(1-r^2)
#signif(qt(qtiles,df=n-2),2)
cairo(file='../../figures/tr.pdf',width=6,height=3)
set.seed(1)
df <- n-2
pars(mfrow=c(1,2))
x <- seq(from=qt(0.001,df=df),to=qt(0.999999,df=df),length.out=100)
y <- dt(x,df=df)
plot(x,y,type='l',bty='n',xlab='t',ylab='f(t)')
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,1))
lines(rep(tt,2),range(y),lty=2)
ltail <- (x<qt(0.025,df=df))
polygon(c(x[ltail][1],x[ltail],tail(x[ltail],1)),c(0,y[ltail],0),col='black')
utail <- (x>qt(0.975,df=df))
polygon(c(x[utail][1],x[utail],tail(x[utail],1)),c(0,y[utail],0),col='black')
Y <- pt(x,df=df)
plot(x,Y,type='l',bty='n',xlab='t',ylab='F(t)')
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(1.5,1))
lines(rep(tt,2),range(Y),lty=2)
lines(range(x),rep(pt(tt,df=df),2),lty=2)
lines(range(x),rep(0.025,2),lty=3)
lines(range(x),rep(0.975,2),lty=3)
lines(rep(-x[utail][1],2),c(0,1),lty=3)
lines(rep(x[utail][1],2),c(0,1),lty=3)
dev.off()

cairo(file='../../figures/RbSr.pdf',width=3,height=3)
pars()
xlim <- c(0,max(RbSr['RbSr',]+2*RbSr['errRbSr',]))
ylim <- c(0.65,max(RbSr['SrSr',]+2*RbSr['errSrSr',]))
plot(RbSr['RbSr',],RbSr['SrSr',],pch=21,bg='white',
     xlim=xlim,ylim=ylim,bty='n',axes=FALSE,
     xlab=expression(''^87*'Rb/'^86*'Sr'),ylab='')
axis(side=1,at=c(0,2,4,6,8,9))
axis(side=2,at=c(0.65,0.7,0.75,0.8,0.85))
mtext(side=2,text=expression(''^87*'Sr/'^86*'Sr'),line=1.2)
dev.off()

tryRbSr <- function(RbSr,a,b,ylim=range(RbSr['SrSr',]),...){
    plot(RbSr['RbSr',],RbSr['SrSr',],type='n',
         xlim=xlim,ylim=ylim,bty='n',axes=FALSE,
         xlab=expression(''^87*'Rb/'^86*'Sr'),ylab='')
    axis(side=1,at=c(0,2,4,6,8,9))
    axis(side=2,at=c(0.6,0.65,0.7,0.75,0.8,0.85,0.9))
    mtext(side=2,text=expression(''^87*'Sr/'^86*'Sr'),line=1.2,cex=0.7)
    lines(xlim,a+b*xlim)
    for (i in 1:n){
        lines(rep(RbSr['RbSr',i],2),c(RbSr['SrSr',i],a+b*RbSr['RbSr',i]),lty=2)
    }
    points(RbSr['RbSr',],RbSr['SrSr',],pch=21,bg='white')
    ss <- sum((a+b*RbSr['RbSr',]-RbSr['SrSr',])^2)
    legend("bottomright",
           c(as.expression(bquote(beta[0] == .(signif(a,3)))),
             bquote(beta[1] == .(signif(b,3))),
             bquote(ss == .(signif(ss,2)))),cex=1.2,bty='n')
}

cairo(file='../../figures/tryRbSr.pdf',width=8,height=2.5)
pars(mar=c(2.7,2.6,0.2,0),mfrow=c(1,3))
tryRbSr(RbSr=RbSr,a=0.6,b=0.03,ylim=c(0.6,0.9))
legend('topleft','a)',bty='n',cex=1.2,adj=c(1,0))
tryRbSr(RbSr=RbSr,a=0.77,b=0,ylim=c(0.6,0.9))
legend('topleft','b)',bty='n',cex=1.2,adj=c(1,0))
tryRbSr(RbSr=RbSr,a=0.698,b=0.0138,ylim=c(0.6,0.9))
legend('topleft','c)',bty='n',cex=1.2,adj=c(1,0))
dev.off()

envelope <- function(X,Y,xlab,ylab,xlim,ylim,predict=FALSE,...){
    fit <- lm(Y ~ X)
    b0 <- fit$coef[1]
    b1 <- fit$coef[2]
    n <- length(X)
    dy2 <- (b0 + b1*X - Y)^2
    vy <- sum(dy2)/(n-2)
    sy <- sqrt(vy)
    den <- sum(X^2-mean(X)^2)
    Eb <- (sy^2)*rbind(c(sum(X^2)/n,-mean(X)),c(-mean(X),1))/den
    tb <- qt(0.025,n-2)
    nn <- 20
    plot(X,Y,type='n',bty='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
    x <- seq(from=xlim[1],to=xlim[2],length.out=nn)
    y <- 0*x
    yl <- 0*x
    yu <- 0*x
    ylp <- 0*x
    yup <- 0*x
    for (i in 1:nn){
        J <- matrix(c(1,x[i]),1,2)
        syp <- sqrt(Eb[1,1] + Eb[2,2]*x[i]^2 + 2*x[i]*Eb[1,2])
        y[i] <- b0 + b1*x[i]
        yl[i] <- y[i] - tb*syp
        yu[i] <- y[i] + tb*syp
        ylp[i] <- y[i] - tb*sqrt(sy^2+syp^2)
        yup[i] <- y[i] + tb*sqrt(sy^2+syp^2)
    }
    polygon(x=c(x,rev(x)),y=c(yu,rev(yl)),col='gray70',border=NA)
    lines(x,y)
    if (predict){
        lines(x,yup,lty=2)
        lines(x,ylp,lty=2)   
    }
    points(X,Y,...)
}

cairo(file='../../figures/envelope.pdf',width=3,height=3)
RbSr <- RbSrGenerator(n=8)
X <- RbSr['RbSr',]
Y <- RbSr['SrSr',]
pars(mar=c(2.5,2.8,0.0,0.0))
envelope(X=X,Y=Y,xlab=expression('x = '^87*'Rb/'^86*'Sr'),
         ylab=expression('y = '^87*'Sr/'^86*'Sr'),
         xlim=c(0,10),ylim=c(0.65,0.87),pch=21,bg='white',cex=0.8)
dev.off()

cairo(file='../../figures/prediction-interval.pdf',width=6,height=2)
pars(mfrow=c(1,3))
RbSr <- RbSrGenerator(n=8)
envelope(X=RbSr['RbSr',],Y=RbSr['SrSr',],
         xlab='x',ylab='y',xlim=c(0,10),
         ylim=c(0.65,0.87),predict=TRUE,
         pch='.',col='gray30')
legend('topleft',legend='a) n = 8',bty='n')
RbSr <- RbSrGenerator(n=100)
envelope(X=RbSr['RbSr',],Y=RbSr['SrSr',],
         xlab='x',ylab='y',xlim=c(0,10),
         ylim=c(0.65,0.87),predict=TRUE,
         pch='.',col='gray30')
legend('topleft',legend='b) n = 100',bty='n')
RbSr <- RbSrGenerator(n=1000)
envelope(X=RbSr['RbSr',],Y=RbSr['SrSr',],
         xlab='x',ylab='y',xlim=c(0,10),
         ylim=c(0.65,0.87),predict=TRUE,
         pch='.',col='gray30')
legend('topleft',legend='c) n = 1000',bty='n')
dev.off()

nr <- 3
nc <- 6
cairo(file='../../figures/regression-p-hacking.pdf',width=10,height=10*nr/nc)
m <- rbind(
    c(1,rep(2,nc),3),
    c(1,4+(1:nc),3),
    c(1,4+((nc+1):(2*nc)),3),
    c(1,4+((2*nc+1):(3*nc)),3),
    c(1,rep(4,nc),3)
)
pars(mar=rep(0,4),mgp=c(1.4,0.5,0))
layout(m,widths=c(0.03,rep(0.97/nc,nc),0.001),
       heights=c(0.02,rep(0.9/nr,nr),0.08))#, layout.show(n=nr*nc+3)
set.seed(5)
n <- 10
plot.new()
plot.new()
plot.new()
plot.new()
for (i in 1:nr){
    for (j in 1:nc){
        X <- runif(10)
        Y <- runif(10)
        plot(X,Y,xlim=c(0,1.0),ylim=c(0,1.3),asp=1,xaxt='n',yaxt='n',type='n')
        fit <- lm(Y ~ X)
        smry <- summary(fit)
        r2 <- signif(smry$r.squared,2)
        pval <- signif(smry$coefficients[2,4],2)
        mtext(bquote('r'^2*'='*.(r2)*', p='*.(pval)),line=-1.6,cex=0.8)
        if (pval < 0.05){
            lines(c(0,1),predict(fit,newdata=data.frame(X=c(0,1))))
        }
        points(X,Y,pch=21,bg='white')
        if (i==nr){
            axis(side=1,at=c(0,0.5,1))
            mtext('x',side=1,line=1.5)
        }
        if (j==1){
            axis(side=2,at=c(0,0.5,1))
            mtext('y',side=2,line=1.5)
        }
    }
}
dev.off()

cairo(file='../../figures/regression-outlier.pdf',width=3,height=3)
pars()
set.seed(5)
x <- c(0.5+runif(10),5)
y <- c(0.5+runif(10),5)
xlim <- c(0,5)
plot(x,y,bty='n',type='n',xlim=xlim,ylim=c(0,5))
fit <- lm(y ~ x)
smry <- summary(fit)
r2 <- signif(smry$r.squared,2)
pval <- signif(smry$coefficients[2,4],2)
lines(xlim,fit$coef[1]+fit$coef[2]*xlim)
points(x,y,pch=21,bg='black')
mtext(bquote('r'^2*'='*.(r2)*', p='*.(pval)),line=-1.6)
dev.off()

cairo(file='../../figures/spurious.pdf',width=9,height=1.8)
pars(mfrow=c(1,5),mar=c(2.3,2.2,1.6,0.2))
set.seed(2)
n <- 50
x <- 100+rnorm(n)
y <- 100+rnorm(n)
z <- 100+rnorm(n)*10
plot(x,y,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('a) r'^2*'='*.(signif(cor(x,y),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('x',side=1,line=1.3,cex=1.0)
mtext('y',side=2,line=1.5,cex=1.0)
plot(x,z,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('b) r'^2*'='*.(signif(cor(x,z),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('x',side=1,line=1.3,cex=1.0)
mtext('z',side=2,line=1.4,cex=1.0)
plot(y,z,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('c) r'^2*'='*.(signif(cor(y,z),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('y',side=1,line=1.3,cex=1.0)
mtext('z',side=2,line=1.4,cex=1.0)
plot(x/z,y/z,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('d) r'^2*'='*.(signif(cor(x/z,y/z),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('x/z',side=1,line=1.3,cex=1.0)
mtext('y/z',side=2,line=1.4,cex=1.0)
plot(x/z,z,pch=21,bg='white',axes=FALSE,xlab='',ylab='')
mtext(bquote('e) r'^2*'='*.(signif(cor(x/z,z),2))),line=0,adj=0,cex=0.9)
axis(1)
axis(2)
mtext('x/z',side=1,line=1.3,cex=1.0)
mtext('z',side=2,line=1.4,cex=1.0)
dev.off()

cairo(file='../../figures/errorfit.pdf',width=3.5,height=3.5)
pars(mar=c(2.5,2.3,0.5,0.3))
X <- c(10,20,30)
Y <- c(20,30,40)
x <- c(10.5,19.5,25.1)
sx <-c(1,1,3)
y <- c(20.5,29.9,48.2)
sy <- c(1,1,5)
rxy <- c(0.9,0.9,-0.9)
tab <- cbind(x,sx,y,sy,rxy)
IsoplotR::scatterplot(tab,ellipse.fill=NA)
mtext('x,X',side=1,line=1.5)
mtext('y,Y',side=2,line=1.5)
yfit <- IsoplotR::york(tab)
fit <- lm(y ~ x)
#abline(a=yfit$a[1],b=yfit$b[1],lty=1)
abline(fit,lty=2)
abline(a=10,b=1,lty=3)
points(x,y,pch=21,bg='black')
points(X,Y,pch=22,bg='white')
text(x,y,labels=1:3,pos=c(3,1,3))
dev.off()

cairo(file='../../figures/yorkfit.pdf',width=7,height=3.5)
pars(mar=c(2.5,2.3,2.0,0.3),mfrow=c(1,2))
X <- c(10,20,30)
Y <- c(20,30,40)
x <- c(10.5,19.5,25.1)
sx <-c(1,1,3)
y <- c(20.5,29.9,45.2)
sy <- c(1,1,5)
rxy <- c(0.9,0.9,-0.9)
tab <- cbind(x,sx,y,sy,rxy)
IsoplotR::scatterplot(tab,ellipse.fill=NA)
mtext('x,X',side=1,line=1.5)
mtext('y,Y',side=2,line=1.5)
yfit <- IsoplotR::york(tab)
fit <- lm(y ~ x)
abline(a=yfit$a[1],b=yfit$b[1],lty=1)
abline(fit,lty=2)
abline(a=10,b=1,lty=3)
points(x,y,pch=21,bg='black')
points(X,Y,pch=22,bg='white')
text(x,y,labels=1:3,pos=c(3,1,3))
legend('topleft',legend='a)',bty='n',adj=c(2,0),cex=1.2)
#
colnames(tab) <- c('X','sX','Y','sY','rXY')
IsoplotR::scatterplot(tab,ellipse.fill=NA,
                      xlim=c(17,33),ylim=c(32,60))
abline(a=yfit$a[1],b=yfit$b[1],lty=1)
abline(fit,lty=2)
abline(a=10,b=1,lty=3)
points(x,y,pch=21,bg='black')
points(X,Y,pch=22,bg='white')
xy <- IsoplotR:::get.york.xy(tab,a=yfit$a[1],b=yfit$b[1])
points(xy,pch=21,bg='white')
lines(rep(x[3],2),c(y[3],100),lty=2,col='gray50')
lines(rep(xy[3,1],2),c(xy[3,2],100),lty=2,col='gray50')
lines(rep(X[3],2),c(Y[3],100),lty=2,col='gray50')
axis(side=3,at=c(x[3],xy[3,1],X[3]),
     labels=c(expression('x'[3]),
              expression(bold('x')[3]),
              expression('X'[3])))
mtext('x,X',side=1,line=1.5)
mtext('y,Y',side=2,line=1.5)
legend('topleft',legend='b)',bty='n',adj=c(2,0),cex=1.2)
dev.off()

cairo(file='../../figures/recentquakes.pdf',width=4,height=2.5)
pars()
ticks <- c(4.5,5,6,7,8)
quakes <- read.csv('recentquakes.csv',header=TRUE)
hist(quakes$mag,breaks=20,main='',xlab='magnitude',xaxt='n',xpd=NA)
axis(side=1,at=ticks,labels=ticks)
dev.off()

cairo(file='../../figures/recentlogquakes.pdf',width=4,height=2.5)
pars()
hist(log(quakes$mag),breaks=20,main='',
     xlab='ln[magnitude]',xpd=NA,ylim=c(0,12500))
dev.off()

# from https://rspatial.org/raster/cases/2-coastline.html
library(raster)
uk <- raster::getData('GADM', country='GBR', level=0)
prj <- paste0("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 ",
              "+x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m")
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
guk <- spTransform(uk, CRS(prj))
duk <- disaggregate(guk)
a <- raster::area(duk)
i <- which.max(a)
b <- duk[i,]
measure_with_ruler <- function(pols, length, lonlat=FALSE) {
    # some sanity checking
    stopifnot(inherits(pols, 'SpatialPolygons'))
    stopifnot(length(pols) == 1)
    # get the coordinates of the polygon
    g <- geom(pols)[, c('x', 'y')]
    nr <- nrow(g)
    # we start at the first point
    pts <- 1
    newpt <- 1
    while(TRUE) {
        # start here
        p <- newpt
        # order the points
        j <- p:(p+nr-1)
        j[j > nr] <- j[j > nr] - nr
        gg <- g[j,]
        # compute distances
        pd <- pointDistance(gg[1,], gg, lonlat)
        # get the first point that is past the end of the ruler
        # this is precise enough for our high resolution coastline
        i <- which(pd > length)[1]
        if (is.na(i)) {
            stop('Ruler is longer than the maximum distance found')
        }
        # get the record number for new point in the original order
        newpt <- i + p
        # stop if past the last point
        if (newpt >= nr) break
        pts <- c(pts, newpt)
    }
    # add the last (incomplete) stick.
    pts <- c(pts, 1)
    # return the locations
    g[pts, ]
}
y <- list()
rulers <- c(10,20,50,100,200) # km
for (i in 1:length(rulers)) {
    y[[i]] <- measure_with_ruler(b, rulers[i]*1000)
}
nrulers <- sapply(y, nrow)
L <- nrulers * rulers

cairo(file='../../figures/britain.pdf',width=7,height=2)
pars(mfrow=c(1,5),mar=c(1,0,0,0))
lat <- 6e5
lon <- 5e5
for (i in length(rulers):1){
    plot(y[[i]],type='l',asp=1,bty='n',xaxt='n',yaxt='n',ann=FALSE)
    arrows(lon,lat,lon,lat+rulers[i]*1000,code=3,angle=90,length=.03,lwd=1)
    text(lon,lat+rulers[i]*500,labels=paste0(rulers[i],' km'),pos=4,xpd=NA)
    mtext(side=1,text=paste0('length = ',L[i],'km'),line=0,cex=0.7)
}
dev.off()

y <- list()
logrulers <- seq(from=log(5),to=log(150),length.out=10)
rulers <- exp(logrulers)
for (i in 1:length(rulers)) {
    y[[i]] <- measure_with_ruler(b, rulers[i]*1000)
}
nrulers <- sapply(y, nrow)
L <- nrulers * rulers

cairo(file='../../figures/loglogbritain.pdf',width=3,height=3)
pars()
fit <- lm(log(L)~log(rulers))
plot(rulers,L,log='xy',bty='n',type='n',
     xlab='length of measuring rod [km]',ylab='length of coast [km]')
lines(rulers,exp(fit$coef[1]+fit$coef[2]*log(rulers)))
points(rulers,L,pch=21,bg='white')
legend('topright',paste0('ln[y] = ',signif(fit$coef[1],3),
       signif(fit$coef[2],3),' ln[x]'),bty='n')
dev.off()

cairo(file='../../figures/fractaldimbritain.pdf',width=3,height=3)
pars()
fit <- lm(log(nrulers)~log(rulers))
plot(rulers,nrulers,log='xy',bty='n',type='n',
     xlab='length of segment [km]',ylab='number of segments')
lines(rulers,exp(fit$coef[1]+fit$coef[2]*log(rulers)))
points(rulers,nrulers,pch=21,bg='white')
legend('topright',paste0('ln[y] = ',signif(fit$coef[1],3),
       signif(fit$coef[2],4),' ln[x]'),bty='n')
dev.off()

cairo(file='../../figures/gutenberg.pdf',width=3,height=3)
pars(mgp=c(1.2,0.5,0))
geostats::gutenberg(quakes$mag,pch=21,bg='white')
dev.off()

cairo(file='../../figures/Finland.pdf',width=3,height=3)
pars()
Finland <- read.csv('Finland.csv',header=TRUE)
n <- 10
sizes <- exp(seq(from=log(2*min(Finland$Lake_area)),
                 to=log(max(Finland$Lake_area)/8),length.out=n))
counts <- rep(0,n)
for (i in 1:n){
    counts[i] <- sum(Finland$Lake_area>=sizes[i])
}
plot(x=sizes,y=counts,log='xy',type='n',
     bty='n',xlab=expression('area (km'^2*')'),
     ylab='',xaxt='n',yaxt='n')
axis(side=1,at=c(0.2,1,10,100,1000),labels=c(0.2,1,10,100,1000))
axis(side=2,at=c(1,10,100,1000),labels=c(1,10,100,1000))
fit <- lm(log(counts) ~ log(sizes))
lines(range(sizes),exp(fit$coef[1]+fit$coef[2]*log(range(sizes))))
points(x=sizes,y=counts,pch=21,bg='white')
mtext('# lakes',side=2,line=1.2)
legend('topright',paste0('y = ',signif(fit$coef[1],3),
                         signif(fit$coef[2],3),' x'),bty='n')
dev.off()

if (FALSE){
    pars(mar=c(3,0,0,0))
    h <- 0.1
    pol <- read.csv('polarity2.csv',header=TRUE)
    nr <- nrow(pol)
    plot(c(pol[1,'start'],pol[nr,'end']),c(0,h),type='n',
         bty='n',axes=FALSE,xlab='time (Ma)',ylab='')
    axis(side=1)
    for (i in 1:nr){
        polygon(x=c(pol[i,'start'],pol[i,'start'],pol[i,'end'],
                    pol[i,'end'],pol[i,'start']),
                y=c(0,h,h,0,0),col='black')
    }
    l <- pol[1:nr,'end']-pol[1:nr,'start']
    durations <- exp(seq(from=log(min(l)),to=log(max(l)),length.out=n))
    counts <- rep(0,n)
    for (i in 1:n){
        counts[i] <- sum(l>=durations[i])
    }
    plot(x=durations,y=counts,log='xy',#type='n',
         bty='n',xlab=expression('length (Ma)'),
         ylab='',xaxt='n',yaxt='n')
}

raster2dat <- function(fname){
    library(raster)
    dat <- raster(fname)
    mat <- t(apply(as.matrix(dat), 2, rev))
    mat[mat<128] <- 1
    mat[mat>=128] <- 0
    mat
}
count_boxes <- function(mat,boxside,nticks){
    mat4plot <- mat
    nboxes <- nticks^2
    for(j in 1:nticks){
        row_id <- ((j-1)*boxside+1):(j*boxside)
        for(k in 1:nticks){
            col_id <- ((k-1)*boxside+1):(k*boxside)
            if(sum(mat[row_id,col_id]) == 0){
                nboxes <- nboxes - 1
            } else {
                mat4plot[row_id,col_id] <- 1
            }
        }
    }
    list(nboxes=nboxes,mat=mat4plot)
}
plotboxmap <- function(mat,boxside,nticks,i,
                       crop=rep(0,4),col=c('white','black'),...){
    count <- count_boxes(mat,boxside=boxside[i],nticks=nticks[i])
    m <- count$mat
    rows2keep <- 1:nrow(m)
    cols2keep <- 1:ncol(m)
    if (crop[1]>0){
        br <- (nrow(m)-floor(nrow(m)*crop[1])):nrow(m) # bottom row
        rows2keep <- rows2keep[-which(rows2keep%in%br)]
    }
    if (crop[2]>0){
        lc <- 1:(1+floor(ncol(m)*crop[2])) # left column
        cols2keep <- cols2keep[-which(cols2keep%in%lc)]
    }
    if (crop[3]){
        tr <- 1:(1+floor(nrow(m)*crop[3])) # top row
        rows2keep <- rows2keep[-which(rows2keep%in%tr)]
    }
    if (crop[4]){
        rc <- (ncol(m)-floor(ncol(m)*crop[4])):ncol(m) # right column
        cols2keep <- cols2keep[-which(cols2keep%in%rc)]
    }
    m4plot <- m[cols2keep,rows2keep]
    image(m4plot,xaxt='n',yaxt='n',ann=FALSE,col=col,
          asp=(1-sum(crop[c(1,3)]))/(1-sum(crop[c(2,4)])),...)
    invisible(count)
}
fractaldim <- function(mat,nticks){
    nboxes <- nticks^2
    for(i in 1:length(nticks)){
        count <- count_boxes(mat=mat,boxside=boxside[i],nticks=nticks[i])
        nboxes[i] <- count$nboxes
    }
    fit <- lm(log(nboxes) ~ log(boxside))
    plot(log(nboxes) ~ log(boxside),type='n',mgp=c(1.2,0.5,0),
         bty='n',xlab='ln[size of boxes]',
         ylab=expression('ln[number of boxes]'))
    abline(fit)
    points(log(nboxes) ~ log(boxside),pch=21,bg='white')
    legend('topright',bty='n',
           legend=paste0('y = ',signif(fit$coef[1],3),
                         signif(fit$coef[2],3),' x'))
}

tiff(file='Britain.tiff',width=512,height=512)
pars(mar=rep(0,4))
plot(b)
dev.off()

nticks <- rep(2,10)^(0:9)
boxside <- rev(nticks)

Britain <- raster2dat('Britain.tiff')
png(file='../../figures/Britainboxes.png',type='cairo',
    family="serif",pointsize=25,
    width=11.6,height=3,res=300,units='in')
pars(mfrow=c(1,4),mar=rep(0,4))
mat <- Britain
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=5,frame=FALSE)
legend('topright',bty='n',legend=paste0(count$nboxes,'/',nticks[5]^2))
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=6,frame=FALSE)
legend('topright',bty='n',legend=paste0(count$nboxes,'/',nticks[6]^2))
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=8,frame=FALSE)
legend('topright',bty='n',legend=paste0(count$nboxes,'/',nticks[8]^2))
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=10,frame=FALSE)
legend('topright',bty='n',legend=paste0(count$nboxes,'/',nticks[10]^2))
dev.off()

cairo(file='../../figures/Britainboxcounts.pdf',width=3,height=3)
pars()
fractaldim(mat=Britain,nticks=nticks)
dev.off()

Corsica <- raster2dat('Corsica.tif')
png(file='../../figures/Corsica.png',type='cairo',
    family="serif",pointsize=25,width=11,height=3,res=300,units='in')
pars(mfrow=c(1,4),mar=rep(0,4))
mat <- Corsica
crop <- c(0.1,0.1,0.1,0.1)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=7,frame=FALSE,crop=crop)
legend('topleft',bty='n',legend=paste0(count$nboxes,'/',nticks[7]^2))
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=8,frame=FALSE,crop=crop)
legend('topleft',bty='n',legend=paste0(count$nboxes,'/',nticks[8]^2))
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=9,frame=FALSE,crop=crop)
legend('topleft',bty='n',legend=paste0(count$nboxes,'/',nticks[9]^2))
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=10,frame=FALSE,crop=crop)
legend('topleft',bty='n',legend=paste0(count$nboxes,'/',nticks[10]^2))
dev.off()

cairo(file='../../figures/Corsicaboxcounts.pdf',width=3,height=3)
pars()
fractaldim(mat=Corsica,nticks=nticks)
dev.off()

drawSnow <- function(lev, x1, y1, x5, y5){
    if (lev == 0){
        lines(x=c(x1,x5),y=c(y1,y5))
    } else {
        deltaX = x5 - x1
        deltaY = y5 - y1
        x2 = x1 + deltaX / 3
        y2 = y1 + deltaY / 3
        x3 = (0.5 * (x1+x5) + sqrt(3) * (y1-y5)/6)
        y3 = (0.5 * (y1+y5) + sqrt(3) * (x5-x1)/6)
        x4 = x1 + 2 * deltaX /3
        y4 = y1 + 2 * deltaY /3
        drawSnow(lev-1, x1, y1, x2, y2)
        drawSnow(lev-1, x2, y2, x3, y3)
        drawSnow(lev-1, x3, y3, x4, y4)
        drawSnow(lev-1, x4, y4, x5, y5)
    }
}

xlim <- c(0,100)
ylim <- c(0,30)

cairo(file='../../figures/koch1.pdf',width=6.5,height=2)
pars(mar=rep(0,4))
plot(xlim,ylim,type='n',asp=1,bty='n',ann=FALSE,xaxt='n',yaxt='n')
level <- 1
drawSnow(lev=level,0,0,100,0)
dev.off()

cairo(file='../../figures/koch2.pdf',width=6.5,height=2)
pars(mar=rep(0,4))
plot(xlim,ylim,type='n',asp=1,bty='n',ann=FALSE,xaxt='n',yaxt='n')
level <- 2
drawSnow(lev=level,0,0,100,0)
dev.off()

cairo(file='../../figures/koch3.pdf',width=6.5,height=2)
pars(mar=rep(0,4))
plot(xlim,ylim,type='n',asp=1,bty='n',ann=FALSE,xaxt='n',yaxt='n')
level <- 3
drawSnow(lev=level,0,0,100,0)
dev.off()

cairo(file='../../figures/koch6.pdf',width=6.5,height=2)
pars(mar=rep(0,4))
plot(xlim,ylim,type='n',asp=1,bty='n',ann=FALSE,xaxt='n',yaxt='n')
level <- 6
drawSnow(lev=level,0,0,100,0)
dev.off()

tiff(file='koch6.tif',width=512,height=512,units='px')
pars(mar=rep(0,4))
plot(xlim,ylim,type='n',asp=1,bty='n',ann=FALSE,xaxt='n',yaxt='n')
level <- 6
drawSnow(lev=level,0,0,100,0)
dev.off()

kochimg <- raster2dat('koch6.tif')
png(file='../../figures/koch.png',type='cairo',width=12,height=1.3,res=300,units='in')
pars(mfrow=c(1,4),mar=rep(0,4))
mat <- kochimg
crop <- c(0.3,0,0.3,0)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=5,frame=FALSE,crop=crop)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=6,frame=FALSE,crop=crop)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=8,frame=FALSE,crop=crop)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,i=10,frame=FALSE,crop=crop)
dev.off()

cairo(file='../../figures/kochboxcounts.pdf',width=3,height=3)
pars()
fractaldim(mat=kochimg,nticks=rep(2,10)^(0:9))
dev.off()

cairo(file='../../figures/sierpinski.pdf',width=12,height=2.5)
pars(mfrow=c(1,4),mar=rep(1,4))
g <- sierpinski(level=1)
image(g,col=c('white','black'),axes=FALSE,asp=1)
g <- sierpinski(level=2)
image(g,col=c('white','black'),axes=FALSE,asp=1)
g <- sierpinski(level=3)
image(g,col=c('white','black'),axes=FALSE,asp=1)
g <- sierpinski(level=4)
image(g,col=c('white','black'),axes=FALSE,asp=1)
dev.off()

tiff(file='sierpinski.tif',width=512,height=512,units='px')
pars(mar=rep(0,4))
g <- sierpinski(level=4)
image(g,col=c('white','black'),axes=FALSE,asp=1)
dev.off()

sierpimg <- raster2dat('sierpinski.tif')
png(file='../../figures/sierpinski.png',type='cairo',width=15,height=3,res=300,units='in')
pars(mfrow=c(1,4),mar=rep(0,4))
mat <- sierpimg
crop <- rep(0,4)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,
                    i=6,frame=FALSE,crop=crop)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,
                    i=7,frame=FALSE,crop=crop)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,
                    i=8,frame=FALSE,crop=crop)
count <- plotboxmap(mat,boxside=boxside,nticks=nticks,
                    i=10,frame=FALSE,crop=crop)
dev.off()

cairo(file='../../figures/sierpinskiboxcounts.pdf',width=3,height=3)
pars()
fractaldim(mat=sierpimg,nticks=rep(2,10)^(0:9))
dev.off()

cantor <- function(n,Y=0){
    if (n <= 0){
        x <- c(0,1)
        y <- rep(Y,2)
    } else {
        xy <- cantor(n-1,Y=Y)
        x <- (1/3)*c(xy$x,xy$x+2)
        y <- rep(xy$y,2)
    }
    list(x=x,y=y)
}
cantorlines <- function(n,Y){
    xy <- cantor(n=n,Y=Y)
    x <- rbind(xy$x[1:(length(xy$x)-1)],xy$x[-1])
    y <- rbind(xy$y[1:(length(xy$y)-1)],xy$y[-1])
    i <- seq(from=1,to=ncol(x),by=2)
    xl <- x[,i]
    yl <- y[,i]
    matlines(xl,yl,lty=1,col='black',lwd=3,xpd=NA)
    text(0,Y,labels=n,pos=2,xpd=NA)
    invisible(xl)
}

cairo(file='../../figures/cantor.pdf',width=4,height=1.5)
pars(mar=c(0,0.5,0.1,0))
plot(c(0,1),y=c(0,1),type='n',bty='n',ann=FALSE,xaxt='n',yaxt='n',xpd=NA)
cantorlines(n=0,Y=1.00)
cantorlines(n=1,Y=0.75)
cantorlines(n=2,Y=0.50)
cantorlines(n=3,Y=0.25)
cantorlines(n=4,Y=0.00)
dev.off()

cairo(file='../../figures/cantorloglog.pdf',width=3,height=3)
pars()
mult <- 1e4
n <- 10
xy <- cantor(n=n,Y=0)
l <- diff(xy$x)[seq(from=2,to=length(xy$x)-1,by=2)]*mult
sizes <- exp(seq(from=log(min(l)),to=log(max(l)),length.out=n))
counts <- rep(0,n)
for (i in 1:n){
    counts[i] <- sum(l>=sizes[i])
}
plot(sizes,counts,log='xy',bty='n',type='n',axes=FALSE,
     xlab='length of segment',ylab='number of segments')
axis(side=1,at=c(1,10,100,1000))
axis(side=2,at=c(1,10,100,1000))
fit <- lm(log(counts) ~ log(sizes))
X <- range(sizes)
Y <- exp(fit$coef[1] + log(X)*fit$coef[2])
lines(X,Y)
points(sizes,counts,pch=21,bg='white')
legend('topright',legend=paste0('ln[y] = ',signif(fit$coef[1],2),
                                signif(fit$coef[2],2),' ln[x]'),bty='n')
dev.off()

src <- rbind(c(0,0),c(.5,sqrt(.75)),c(1,0))
cairo(file='../../figures/pendulum.pdf',width=3.5,height=3.5)
pars(mar=rep(0,4))
p <- persp(x=c(-2,2),y=c(-2,2),matrix(0,2,2),axes=FALSE,
           zlim=c(0,2),phi=30,theta=-30,zlab='z',border='gray50')
ell <- get.ellipse(mu=1.5*src[3,]-0.7,a=0.3,b=0.15,theta=0)
xyell <- trans3d(x=ell[,1],y=ell[,2],z=rep(0,nrow(ell)),pmat=p)
lines(xyell$x,xyell$y)
text(xyell$x[20],xyell$y[20],labels='3',pos=3,offset=0.1)
ell <- get.ellipse(mu=1.5*src[2,]-0.7,a=0.3,b=0.15,theta=0)
xyell <- trans3d(x=ell[,1],y=ell[,2],z=rep(0,nrow(ell)),pmat=p)
lines(xyell$x,xyell$y)
text(xyell$x[20],xyell$y[20],labels='2',pos=3,offset=0.1)
ell <- get.ellipse(mu=1.5*src[1,]-0.7,a=0.3,b=0.15,theta=0)
xyell <- trans3d(x=ell[,1],y=ell[,2],z=rep(0,nrow(ell)),pmat=p)
lines(xyell$x,xyell$y)
text(xyell$x[20],xyell$y[20],labels='1',pos=3,offset=0.1)
xystring <- trans3d(x=c(-0.8,0),y=c(-.25,0),z=c(0.2,2),pmat=p)
lines(xystring$x,xystring$y)
ell <- get.ellipse(mu=c(xystring$x[1],xystring$y[1]),a=0.02,b=0.02,theta=0)
polygon(ell,col='black')
dev.off()

throw <- function(startpos=c(-2,2),startvel=c(0,0),src,plot=TRUE){
    niter <- 10000
    pos <- matrix(0,niter,2)
    vel <- matrix(0,niter,2)
    acc <- matrix(0,niter,2)
    pos[1,] <- startpos
    pos[2,] <- startpos
    vel[1,] <- startvel
    vel[2,] <- startvel
    kf <- 1
    friction <- .2
    m_fHeight <- .2
    src[,1] <- src[,1] - 0.5
    src[,2] <- src[,2] - 1 + 0.5/cos(30*pi/180)
    dt <- .02
    for (i in 3:niter){
        pos[i,1] <- pos[i-1,1] + vel[i-1,1]*dt +
            ((2/3)*acc[i-1,1] - (1/6)*acc[i-2,1]) * dt^2
        pos[i,2] <- pos[i-1,2] + vel[i-1,2]*dt +
            ((2/3)*acc[i-1,2] - (1/6)*acc[i-2,2]) * dt^2
        D <- Inf
        for (j in 1:nrow(src)){
            r <- pos[i,] - src[j,];
            d <- sqrt( r[1]^2 + r[2]^2 + m_fHeight^2 )
            if (d<D){
                best <- j
                D <- d
            }
            if (FALSE){
                acc[i,1] <- acc[i,1] - kf*r[1]
                acc[i,2] <- acc[i,2] - kf*r[2]
            } else {
                acc[i,1] <- acc[i,1] - kf*r[1]/d^3
                acc[i,2] <- acc[i,2] - kf*r[2]/d^3
            }
        }
        acc[i,1] <- acc[i,1] - vel[i-1,1]*friction
        acc[i,2] <- acc[i,2] - vel[i-1,2]*friction
        vel[i,1] <- vel[i-1,1] +
            ( (1/3)*acc[i,1] + (5/6)*acc[i-1,1] - (1/6)*acc[i-2,1] )*dt
        vel[i,2] <- vel[i-1,2] +
            ( (1/3)*acc[i,2] + (5/6)*acc[i-1,2] - (1/6)*acc[i-2,2] )*dt
        if (sum(abs(vel[i,]))<.001){
            break
        }
    }
    if (plot){
        pars(mar=rep(0,4))
        plot(c(-2,2),c(-2,2),type='n',bty='n',ann=FALSE,xaxt='n',yaxt='n')
        lines(pos[1:i,])
        points(src,pch=21,bg='white',cex=2.5,lwd=2)
        text(src,labels=1:3)
    }
    invisible(best)
}

cairo(file='../../figures/3magnets1.pdf',width=3,height=3)
throw(startpos=c(-2,2),startvel=c(0,-2),src=src,plot=TRUE)
dev.off()

cairo(file='../../figures/3magnets2.pdf',width=3,height=3)
throw(startpos=c(-1.9,2),startvel=c(0,-2),src=src,plot=TRUE)
dev.off()

if (FALSE){
    res <- 512
    mat <- matrix(0,res,res)
    x <- seq(from=-2,to=2,length.out=res)
    y <- seq(from=-2,to=2,length.out=res)
    for (i in 1:res){
        print(i)
        for (j in 1:res){
            mat[i,j] <- throw(startpos=c(x[i],y[j]),src=src,plot=FALSE)
        }
    }
}

png(file='../../figures/3magnets.png',type='cairo',
    family="serif",pointsize=25,width=6,height=6,res=300,units='in')
load('3magnets.RData')
pars(mar=rep(0,4))
image(x,y,mat,col=c('white','gray50','black'))
points(src,pch=21,bg='white',cex=2)
text(src,labels=1:3)
dev.off()

cairo(file='../../figures/LLpois.pdf',width=3,height=3)
pars()
LL <- function(lambda,k){
  k * log(lambda) - lambda - sum(1:k)
}
nl <- 100
l <- seq(from=0,to=20,length.out=nl)
L <- rep(0,nl)
for (i in 1:nl){
  L[i] <- LL(lambda=l[i],k=4)
}
plot(l,L,type='l',xlab=expression(lambda),ylab='LL')
dev.off()

cairo(file='../../figures/PCA2Ddata.pdf',width=2.2,height=2.2)
pars()
X <- matrix(c(-1,3,4,7,2,3),nrow=3,ncol=2)
colnames(X) <- c('a','b')
plot(X,type='n',xlim=c(-1.2,4.2),ylim=c(1.8,7.2))
text(X,labels=1:3)
dev.off()

PCA2D <- function(X,fig=1){
    pc <- prcomp(X)
    # calculate the end-points of two lines marking the principal components (PC):
    CL <- matrix(NA,4,2) # initialise the matrix of coordinates
    CL[1:2,] <- matrix(1,2,1) %*% pc$center + diag(pc$sdev) %*% t(pc$rotation)
    CL[3:4,] <- matrix(1,2,1) %*% pc$center - diag(pc$sdev) %*% t(pc$rotation)
    if (fig==1){ # initialise the 1st panel:
        rx <- range(X[,'a'],CL[,1]) # range of x-values
        ry <- range(X[,'b'],CL[,2]) # range of y-values
        plot(rx,ry,type='n',asp=1,xlab='a',ylab='b',
             xlim=c(-1.2,4.8),ylim=c(1.2,7.2))
        text(X,labels=1:3)
        # draw the line marking the 1st PC:
        lines(CL[c(1,3),])
        text(CL[3,1],CL[3,2],labels='PC1',pos=4)
        # draw the line marking the 2nd PC:
        lines(CL[c(2,4),])
        text(CL[2,1],CL[2,2],labels='PC2',pos=4)
        # add the centre point as a white square:
        points(t(pc$center),pch=22,bg='white')
    }
    if (fig==2){ # initialise the 2nd panel:
        plot(range(pc$x),c(1,3.5),type='n',bty='n',
             xaxt='n',yaxt='n',xlab='',ylab='',ylim=c(1,3.5))
        Axis(side=1)
        # plot the 1st PC scores as a 1D configuration:
        lines(pc$x[,'PC1'],rep(2,3))
        points(pc$x[,'PC1'],rep(2,3),pch=21,bg='white')
        text(pc$x[,'PC1'],rep(2,3),labels=1:3,pos=c(1,1,3))
        text(min(pc$x[,'PC1']),2,labels='PC1',pos=2,xpd=NA)
        # plot the 2nd PC scores as a 1D configuration:
        lines(pc$x[,'PC2'],rep(3,3))
        points(pc$x[,'PC2'],rep(3,3),pch=21,bg='white')
        text(pc$x[,'PC2'],rep(3,3),labels=1:3,pos=1)
        text(min(pc$x[,'PC2']),3,labels='PC2',pos=2)
    }
    if (fig==3){ # plot both PCA scores and the loadings in the 3rd panel:
        biplot(pc,col='black')
    }
    if (fig==4){ # plot the weights of the PCs in the 4th panel:
        w <- pc$sdev^2
        names(w) <- colnames(pc$x)
        barplot(w)
    }
}

pc <- prcomp(X)

cairo(file='../../figures/PCA2D1.pdf',width=2.5,height=2.5)
pars()
PCA2D(X=X,fig=1)
dev.off()

cairo(file='../../figures/PCA2D2.pdf',width=2.5,height=1.5)
pars(mar=c(1.5,1.6,0.5,0.25))
PCA2D(X=X,fig=2)
dev.off()

cairo(file='../../figures/PCA2D3.pdf',width=2.5,height=2.5)
pars(mar=c(2.4,2.3,1.5,1.5))
PCA2D(X=X,fig=3)
dev.off()

cairo(file='../../figures/USArrests.pdf',width=5,height=5)
pars(mar=c(2.5,2.5,1.5,1.5))
pc <- prcomp(USArrests, scale=TRUE)
biplot(pc,col=c('gray50','black'),xlim=c(-0.3,0.33))
dev.off()

cairo(file='../../figures/eurodist.pdf',width=4,height=4)
pars()
ED <- eurodist
cities <- attr(ED,'Labels')
iLyon <- which(cities%in%'Lyons')
cities[iLyon] <- 'Lyon'
attr(ED,'Labels') <- cities
conf <- MASS::isoMDS(ED)
plot(c(-2500,2500),c(1900,-1900),type='n',
     ylim=c(1950,-1950),xlab='x (km)',ylab='y (km)')
text(conf$points,labels=labels(ED))
dev.off()

cairo(file='../../figures/Shepard.pdf',width=3,height=3)
pars()
plot(MASS::Shepard(d=ED,x=conf$points),pch=19,cex=0.7,
     xlab='true distance',ylab='fitted distance')
dev.off()

cairo(file='../../figures/DZmds.pdf',width=4,height=3)
pars(mar=c(2.5,2.5,0,0.5))
library(IsoplotR)
set.seed(1)
dz <- provenance::read.distributional('DZages.csv',check.names=FALSE)$x
DZ <- lapply(dz,jitter)
d <- matrix(0,13,13)
for (i in 1:13){
    for (j in 1:13){
        d[i,j] <- ks.test(DZ[[i]],DZ[[j]])$statistic
    }
}
nms <- c(1:10,'L','T','Y')
rownames(d) <- nms
colnames(d) <- nms
m <- isoMDS(d)
plot(m$points,type='n',xlab='x',ylab='y',asp=1)
text(m$points,labels=nms)
dev.off()

#as.dist(round(100*d),upper=TRUE)

cairo(file='../../figures/Iris.pdf',width=6,height=6)
IRIS <- iris
names(IRIS) <- c("Sepal Length","Sepal Width",
                 "Petal Length","Petal Width","Species")
pch <- c(rep(1,50),rep(2,50),rep(3,50))
plot(IRIS[,1:4],cex=0.8,pch=pch,oma=rep(1,4),mar=rep(0,4),mgp=c(1.2,0.5,0))
dev.off()

cairo(file='../../figures/kmeans1.pdf',width=3,height=3)
pars()
xy <- IRIS[,c("Sepal Width","Petal Width")]
plot(unique(xy),pch=19,xlab='x',ylab='y',cex=0.5)
dev.off()

cairo(file='../../figures/kmeans2.pdf',width=3,height=3)
pars()
i <- c(48,100,130)
centers <- xy[i,]
pch <- rep('?',150)
cex <- rep(0.75,150)
bg <- rep('white',150)
plot(xy,type='n',xlab='x',ylab='y')
marker1 <- which(xy[,1]%in%xy[i[1],1] & xy[,2]%in%xy[i[1],2])
marker2 <- which(xy[,1]%in%xy[i[2],1] & xy[,2]%in%xy[i[2],2])
marker3 <- which(xy[,1]%in%xy[i[3],1] & xy[,2]%in%xy[i[3],2])
other <- (1:150)[-c(marker1,marker2,marker3)]
points(xy[marker1[1],],pch=1,cex=1.2,lwd=2)
points(xy[marker2[1],],pch=2,cex=1.2,lwd=2)
points(xy[marker3[1],],pch=3,cex=1.2,lwd=2)
points(unique(xy[other,]),pch='?',cex=0.8,col='gray50')
dev.off()

cairo(file='../../figures/kmeans3.pdf',width=3,height=3)
pars()
d <- as.matrix(dist(rbind(centers,xy)))
group1 <- (d[1,-(1:3)] <= d[2,-(1:3)]) & (d[1,-(1:3)] <= d[3,-(1:3)])
group2 <- (d[2,-(1:3)] < d[1,-(1:3)]) & (d[2,-(1:3)] <= d[3,-(1:3)])
group3 <- (d[3,-(1:3)] < d[1,-(1:3)]) & (d[3,-(1:3)] < d[2,-(1:3)])
plot(xy,type='n',xlab='x',ylab='y')
points(unique(xy[group1,]),pch=1,cex=0.8)
points(unique(xy[group2,]),pch=2,cex=0.8)
points(unique(xy[group3,]),pch=3,cex=0.8)
dev.off()

cairo(file='../../figures/kmeans4.pdf',width=3,height=3)
pars()
fit <- kmeans(xy,centers=centers,iter.max=1)
plot(xy,type='n',xlab='x',ylab='y')
points(unique(xy[group1,]),pch=1,col='gray50',cex=0.8)
points(unique(xy[group2,]),pch=2,col='gray50',cex=0.8)
points(unique(xy[group3,]),pch=3,col='gray50',cex=0.8)
mean1 <- colMeans(xy[group1,])
mean2 <- colMeans(xy[group2,])
mean3 <- colMeans(xy[group3,])
points(x=mean1[1],y=mean1[2],pch=1,cex=1.2,lwd=2)
points(x=mean2[1],y=mean2[2],pch=2,cex=1.2,lwd=2)
points(x=mean3[1],y=mean3[2],pch=3,cex=1.2,lwd=2)
dev.off()

cairo(file='../../figures/kmeans5.pdf',width=3,height=3)
pars()
set.seed(5)
fit <- kmeans(xy,centers=3)
pch <- fit$cluster
pch[fit$cluster==1] <- 1
pch[fit$cluster==2] <- 2
pch[fit$cluster==3] <- 3
plot(xy,xlab='x',ylab='y',pch=pch,cex=0.8,col='gray50')
points(fit$centers,cex=1.2,lwd=2,pch=c(1,2,3))
dev.off()

set.seed(18)
fit <- kmeans(iris[,-5],centers=3,nstart=20)
#table(fit$cluster,iris[,'Species'])

cairo(file='../../figures/hierarchical1.pdf',width=3,height=3)
pars()
set.seed(1)
ab <- rbind(c(6,90),
            c(91,43),
            c(17,95),
            c(20,70),
            c(80,06))
plot(ab,type='n',asp=1,xlim=c(0,100),ylim=c(0,100),xlab='x',ylab='y')
text(ab,labels=1:5,cex=1.4)
dev.off()

d1 <- signif(dist(ab,upper=TRUE),3)

cairo(file='../../figures/hierarchical2.pdf',width=3,height=3)
pars()
plot(ab,type='n',asp=1,xlim=c(0,100),ylim=c(0,100),xlab='x',ylab='y')
lines(ab[c(1,3),],col='gray50')
text(ab,labels=1:5,cex=1.4)
dev.off()

cairo(file='../../figures/tree1.pdf',width=2.5,height=2.5)
pars(mar=rep(0,4))
plot(c(1,5),c(-0.5,1.2),type='n',ann=FALSE,bty='n',axes=FALSE)
text(x=1:2,y=0,labels=c(1,3))
text(x=3:5,y=0,labels=c(2,4,5))
lines(x=c(1.5,1.5),y=c(1,0.5))
lines(x=c(1,2),y=c(0.5,0.5))
lines(x=c(1,1),y=c(0.5,0.1))
lines(x=c(2,2),y=c(0.1,0.5))
dev.off()

d1b <- as.matrix(d1)[c(1,3,2,4,5),c(1,3,2,4,5)]
d2 <- matrix(0,4,4)
d2[2:4,2:4] <- d1b[3:5,3:5]
for (i in 2:4){
    d2[1,i] <- max(d1b[1,i+1],d1b[2,i+1])
    d2[i,1] <- d2[1,i]
}

cairo(file='../../figures/hierarchical3.pdf',width=3,height=3)
pars()
plot(ab,type='n',asp=1,xlim=c(0,100),ylim=c(0,100),xlab='x',ylab='y')
lines(ab[c(1,3),],col='gray50')
lines(ab[c(1,4),],col='gray50')
lines(ab[c(4,3),],col='gray50')
text(ab,labels=1:5,cex=1.4)
dev.off()

cairo(file='../../figures/tree2.pdf',width=2.5,height=2.5)
pars(mar=rep(0,4))
plot(c(1,5),c(-0.5,1.7),type='n',ann=FALSE,bty='n',axes=FALSE)
text(x=1:2,y=0,labels=c(1,3))
text(x=3,y=0,labels=4)
text(x=4:5,y=0,labels=c(2,5))
lines(x=c(1.5,1.5),y=c(1,0.5))
lines(x=c(1,2),y=c(0.5,0.5))
lines(x=c(1,1),y=c(0.5,0.1))
lines(x=c(2,2),y=c(0.1,0.5))
lines(x=c(2.25,2.25),y=c(1,1.5))
lines(x=c(1.5,3),y=c(1,1))
lines(x=c(3,3),y=c(0.1,1))
dev.off()

d2b <- as.matrix(d2)[c(1,3,2,4),c(1,3,2,4)]
d3 <- matrix(0,3,3)
d3[2:3,2:3] <- d1b[3:4,3:4]
for (i in 2:3){
    d3[1,i] <- max(d2b[1,i+1],d1b[2,i+1])
    d3[i,1] <- d3[1,i]
}

cairo(file='../../figures/hierarchical4.pdf',width=3,height=3)
pars()
plot(ab,type='n',asp=1,xlim=c(0,100),ylim=c(0,100),xlab='x',ylab='y')
lines(ab[c(1,3),],col='gray50')
lines(ab[c(1,4),],col='gray50')
lines(ab[c(4,3),],col='gray50')
lines(ab[c(2,5),],col='gray50')
text(ab,labels=1:5,cex=1.4)
dev.off()

cairo(file='../../figures/tree3.pdf',width=2.5,height=2.5)
pars(mar=rep(0,4))
plot(c(1,5),c(-0.5,1.7),type='n',ann=FALSE,bty='n',axes=FALSE)
text(x=1:2,y=0,labels=c(1,3))
text(x=3,y=0,labels=4)
text(x=4:5,y=0,labels=c(2,5))
y <- 0.5
lines(x=c(1.5,1.5),y=c(1,0.5))
lines(x=c(1,2),y=c(y,y))
lines(x=c(1,1),y=c(y,0.1))
lines(x=c(2,2),y=c(0.1,y))
y <- 1
lines(x=c(2.25,2.25),y=c(y,1.5))
lines(x=c(1.5,3),y=c(y,y))
lines(x=c(3,3),y=c(0.1,y))
lines(x=c(4.5,4.5),y=c(y,1.5))
lines(x=c(4,5),y=c(y,y))
lines(x=c(4,4),y=c(y,0.1))
lines(x=c(5,5),y=c(y,0.1))
dev.off()

cairo(file='../../figures/hierarchicaltree.pdf',width=3,height=3)
pars()
tree <- hclust(d1)
plot(tree,main='',sub='',xlab='')
dev.off()

cairo(file='../../figures/irishtree.pdf',width=8,height=1.5)
p <- par(oma=rep(0,4),mar=c(0,2,0,0),mgp=c(1.2,0.5,0))
tree <- hclust(dist(iris[,-5]))
plot(tree,main='',sub='',xlab='',labels=FALSE,ann=TRUE)
par(p)
dev.off()

cairo(file='../../figures/irishtreecut.pdf',width=8,height=1.5)
p <- par(oma=rep(0,4),mar=c(0,2,0,0),mgp=c(1.2,0.5,0))
plot(tree,main='',sub='',xlab='',labels=FALSE,ann=TRUE)
lines(x=c(1,150),y=rep(3.6,2),lty=2)
par(p)
dev.off()

#treecut <- cutree(tree,k=3)
#table(treecut, iris[,'Species'])

LDApredict <- utils::getFromNamespace("predict.lda", "MASS")
QDApredict <- utils::getFromNamespace("predict.qda", "MASS")
DApredict <- function(fit,dat){
    if (class(fit)%in%'lda'){
        out <- LDApredict(fit,newdata=dat)
    } else {
        out <- QDApredict(fit,newdata=dat)
    }
}
construct_DA <- function(uv,AFFINITY,quadratic=FALSE,plot=FALSE,padding=0,...){
    u <- uv[,1]
    v <- uv[,2]
    nn <- 1000
    ugrid <- seq(from=min(u,na.rm=TRUE)-padding,
                 to=max(u,na.rm=TRUE)+padding,length.out=nn)
    vgrid <- seq(from=min(v,na.rm=TRUE)-padding,
                 to=max(v,na.rm=TRUE)+padding,length.out=nn)
    uvgrid <- expand.grid(ugrid,vgrid)
    out <- list()
    if (quadratic){
        out$fit <- MASS::qda(AFFINITY ~ u + v,na.action='na.omit')
        nt <- 500
    } else {
        out$fit <- MASS::lda(AFFINITY ~ u + v,na.action='na.omit')
        nt <- 50
    }
    out$contours <- list()
    pr <- DApredict(fit=out$fit,dat=data.frame(u=uvgrid[,1],v=uvgrid[,2]))
    z <- matrix(as.numeric(pr$class),nrow=nn,ncol=nn)
    if (plot) plot(ugrid,vgrid,type='n',...)
    contours <- grDevices::contourLines(ugrid,vgrid,z,levels=c(1.5,2.5))
    for (i in seq_along(contours)){
        nx <- length(contours[[i]]$x)
        j <- seq(from=1,to=nx,length.out=nt)
        x <- contours[[i]]$x[j]
        y <- contours[[i]]$y[j]
        out$contours[[i]] <- cbind(x,y)
        if (plot) graphics::lines(x,y)
    }
    invisible(out)
}

cairo(file='../../figures/QDA.pdf',width=7.5,height=2.5)
pars(mfrow=c(1,3))
set.seed(2)
n <- 100
m1 <- c(-1.2,2)
m2 <- c(1.5,3)
m3 <- c(0,-2)
s1 <- rbind(c(1,0.8),c(0.8,2))
s2 <- rbind(c(0.8,0.5),c(0.5,1.2))
s3 <- rbind(c(2,-0.9),c(-0.9,1))
xy1 <- MASS::mvrnorm(n=n,mu=m1,Sigma=s1)
xy2 <- MASS::mvrnorm(n=n,mu=m2,Sigma=s2)
xy3 <- MASS::mvrnorm(n=n,mu=m3,Sigma=s3)
xy <- rbind(xy1,xy2,xy3)
x <- xy[,1]
y <- xy[,2]
# panel 1
groups <- c(rep(1,n),rep(2,n),rep(3,n))
plot(x,y,pch=groups,xlab='x',ylab='y',cex=0.7)
legend('topleft',legend='a)',bty='n',adj=c(2,0),cex=1.2)
# panel 2
ell1 <- IsoplotR::ellipse(x=m1[1],y=m1[2],covmat=s1)
ell2 <- IsoplotR::ellipse(x=m2[1],y=m2[2],covmat=s2)
ell3 <- IsoplotR::ellipse(x=m3[1],y=m3[2],covmat=s3)
plot(x,y,pch=groups,xlab='x',ylab='y',cex=0.7,col='gray50')
polygon(ell1,border='black',col=NULL)
polygon(ell2,border='black',col=NULL)
polygon(ell3,border='black',col=NULL)
legend('topleft',legend='b)',bty='n',adj=c(2,0),cex=1.2)
# panel 3
qd <- construct_DA(xy,AFFINITY=groups,quadratic=TRUE,
                   plot=TRUE,xlab='x',ylab='y')
points(x,y,pch=groups,xlab='x',ylab='y',col='gray50',cex=0.7)
legend('topleft',legend='c)',bty='n',adj=c(2,0),cex=1.2)
ell <- get.ellipse(mu=m1,a=0.02,b=0.02,theta=0)
dev.off()

cairo(file='../../figures/LDA.pdf',width=7.5,height=2.5)
pars(mfrow=c(1,3))
set.seed(2)
n <- 100
m1 <- c(-1.2,2)
m2 <- c(3,4)
m3 <- c(0,-2)
s1 <- rbind(c(1,0.8),c(0.8,2))
xy1 <- MASS::mvrnorm(n=n,mu=m1,Sigma=s1)
xy2 <- MASS::mvrnorm(n=n,mu=m2,Sigma=s1)
xy3 <- MASS::mvrnorm(n=n,mu=m3,Sigma=s1)
xy <- rbind(xy1,xy2,xy3)
x <- xy[,1]
y <- xy[,2]
# panel 1
groups <- c(rep(1,n),rep(2,n),rep(3,n))
plot(x,y,pch=groups,xlab='x',ylab='y',cex=0.7)
legend('topleft',legend='a)',bty='n',adj=c(2,0),cex=1.2)
# panel 2
ell1 <- IsoplotR::ellipse(x=m1[1],y=m1[2],covmat=s1)
ell2 <- IsoplotR::ellipse(x=m2[1],y=m2[2],covmat=s1)
ell3 <- IsoplotR::ellipse(x=m3[1],y=m3[2],covmat=s1)
plot(x,y,pch=groups,xlab='x',ylab='y',cex=0.7,col='gray50')
polygon(ell1,border='black',col=NULL)
polygon(ell2,border='black',col=NULL)
polygon(ell3,border='black',col=NULL)
legend('topleft',legend='b)',bty='n',adj=c(2,0),cex=1.2)
# panel 3
qd <- construct_DA(xy,AFFINITY=groups,quadratic=FALSE,
                   plot=TRUE,xlab='x',ylab='y')
points(x,y,pch=groups,xlab='x',ylab='y',col='gray50',cex=0.7)
legend('topleft',legend='c)',bty='n',adj=c(2,0),cex=1.2)
ell <- get.ellipse(mu=m1,a=0.02,b=0.02,theta=0)
dev.off()

cairo(file='../../figures/PCAvsLDA.pdf',width=7.5,height=2.5)
pars(mfrow=c(1,3))
set.seed(3)
n <- 100
m1 <- c(-0.5,0.5)
m2 <- c(0.5,-0.5)
s <- rbind(c(1,0.95),c(0.95,1))
groups <- c(rep(1,n),rep(3,n))
xy1 <- MASS::mvrnorm(n=n,mu=m1,Sigma=s)
xy2 <- MASS::mvrnorm(n=n,mu=m2,Sigma=s)
xy <- rbind(xy1,xy2)
# panel 1
plot(xy,pch=groups,asp=1,xlab='x',ylab='y',ylim=c(-3,3))
legend('topleft',legend='a)',bty='n',adj=c(2,0),cex=1.2)
# panel 2
plot(xy,pch=groups,asp=1,col='gray50',xlab='x',ylab='y',ylim=c(-3,3))
legend('topleft',legend='b)',bty='n',adj=c(2,0),cex=1.2)
mu <- colMeans(xy)
Sigma <- cov(xy)
ell <- IsoplotR::ellipse(x=mu[1],y=mu[2],covmat=Sigma)
polygon(ell,border='black',col=NULL)
pc <- princomp(xy)
lines(rbind(10*pc$loadings[,1],-10*pc$loadings[,1]))
# panel 3
plot(xy,pch=groups,asp=1,col='gray50',xlab='x',ylab='y',ylim=c(-3,3))
legend('topleft',legend='c)',bty='n',adj=c(2,1),cex=1.2)
ell1 <- IsoplotR::ellipse(x=m1[1],y=m1[2],covmat=s)
polygon(ell1,border='black',col=NULL)
ell2 <- IsoplotR::ellipse(x=m2[1],y=m2[2],covmat=s)
polygon(ell2,border='black',col=NULL)
ld <- lda(groups ~ xy[,1] + xy[,2])
lines(cbind(unlist(2*ld$scaling),-unlist(2*ld$scaling)))
dev.off()

cairo(file='../../figures/LDAiris.pdf',width=4,height=4)
pars()
ldiris <- lda(Species ~ ., data=iris)
cldiris <- coef(ldiris)
pch <- c(rep(1,50),rep(2,50),rep(3,50))
pred <- predict(ldiris)
ld <- construct_DA(pred$x,AFFINITY=iris[,'Species'],
                   quadratic=FALSE,plot=TRUE,
                   xlab='LD1',ylab='LD2')
text(-3,2.7,'virginica',pos=2,col='gray50',srt=85)
text(0.4,2.7,'versicolor',pos=2,col='gray50',srt=90)
text(4,2.7,'setosa',pos=2,col='gray50',srt=95)
points(pred$x,pch=pch)
dev.off()

cairo(file='../../figures/LDAnewiris.pdf',width=4,height=4)
pars()
pred <- predict(ldiris)
ld <- construct_DA(pred$x,AFFINITY=iris[,'Species'],
                   quadratic=FALSE,plot=TRUE,
                   xlab='LD1',ylab='LD2')
text(-3,2.7,'virginica',pos=2,col='gray50',srt=85)
text(0.4,2.7,'versicolor',pos=2,col='gray50',srt=90)
text(4,2.7,'setosa',pos=2,col='gray50',srt=95)
newflower <- data.frame(Sepal.Length=6.0,Sepal.Width=3.0,
                        Petal.Length=5.0,Petal.Width=1.5)
pred <- predict(ldiris,newdata=newflower)
points(pred$x[1],pred$x[2],pch=8)
dev.off()

means <- colSums(ldiris$prior*ldiris$means)
scaling <- ldiris$scaling
x <- scale(newflower, center = means, scale = FALSE) %*% scaling

library(rpart)
CARTpredict <- utils::getFromNamespace("predict.rpart", "rpart")
CART <- function(u,v,tree,lwd=1,...){
    nn <- 2000
    padding <- 0
    ugrid <- seq(from=min(u,na.rm=TRUE)-padding,
                 to=max(u,na.rm=TRUE)+padding,length.out=nn)
    vgrid <- seq(from=min(v,na.rm=TRUE)-padding,
                 to=max(v,na.rm=TRUE)+padding,length.out=nn)
    uvgrid <- expand.grid(ugrid,vgrid)
    nt <- 500
    pr <- CARTpredict(object=tree,newdata=data.frame(x=uvgrid[,1],y=uvgrid[,2]))
    klasse <- apply(pr,1,function(x){which(x>0.5)})
    z <- matrix(klasse,nrow=nn,ncol=nn)
    plot(ugrid,vgrid,type='n',...)
    contouren <- grDevices::contourLines(ugrid,vgrid,z,levels=c(1.5,2.5))
    for (i in seq_along(contouren)){
        nx <- length(contouren[[i]]$x)
        j <- seq(from=1,to=nx,length.out=nt)
        x <- contouren[[i]]$x[j]
        y <- contouren[[i]]$y[j]
        graphics::lines(x,y,lwd=lwd)
    }
}

cairo(file='../../figures/CARTdata.pdf',width=3.5,height=3.5)
pars()
set.seed(4)
lwd <- 1
n <- 150
m1 <- c(1.5,3)
m2 <- c(-3,0)
m3 <- c(1.5,-3)
m4 <- c(1.5,0)
m5 <- c(-6,0)
mult <- 1
s1 <- rbind(c(4,0),c(0,1))*mult
s2 <- rbind(c(1,0),c(0,4))*mult
s3 <- rbind(c(4,0),c(0,1))*mult
s4 <- rbind(c(2,0),c(0,1))*mult
s5 <- rbind(c(1,0),c(0,4))*mult
xy1 <- MASS::mvrnorm(n=n/3,mu=m1,Sigma=s1)
xy2 <- MASS::mvrnorm(n=n/3,mu=m2,Sigma=s2)
xy3 <- MASS::mvrnorm(n=n/3,mu=m3,Sigma=s3)
xy4 <- MASS::mvrnorm(n=n/2,mu=m4,Sigma=s4)
xy5 <- MASS::mvrnorm(n=n/2,mu=m5,Sigma=s5)
xy <- rbind(xy1,xy2,xy3,xy4,xy5)
x <- xy[,1]
y <- xy[,2]
mx <- min(x)
Mx <- max(x)
my <- min(y)
My <- max(y)
dx <- diff(range(x))
dy <- diff(range(y))
groups <- c(rep(1,n),rep(2,n))
col <- c(rep('white',n),rep('black',n))
plot(x,y,pch=21,axes=TRUE,ann=TRUE,lwd=lwd,xlab='x',ylab='y',bg=col)
dev.off()

plotframe <- function(x,y){
    mx <- min(x)
    Mx <- max(x)
    my <- min(y)
    My <- max(y)
    lines(x=range(x),y=rep(my,2))
    lines(x=range(x),y=rep(My,2))
    lines(x=rep(mx,2),y=range(y))
    lines(x=rep(Mx,2),y=range(y))
}

cairo(file='../../figures/Qtries.pdf',width=10.5,height=3.5)
pars(mfrow=c(1,3),mar=c(2.5,2.3,1.0,0))
col <- c(rep('white',n),rep('gray50',n))
plot(x,y,pch=21,axes=TRUE,ann=TRUE,col='gray50',
     bty='n',lwd=lwd,xlim=c(mx-dx/50,Mx+dx/50),
     ylim=c(my-dy/50,My+dy/50),bg=col,
     xaxs="i",yaxs="i",xlab='x',ylab='y')
plotframe(x,y)
s <- 1
lines(x=rep(s,2),y=range(y))
p1 <- sum(x[1:n] < s)/n
p2 <- sum(x[(n+1):(2*n)] < s)/n
Q <- signif(p1*(1-p1) + p2*(1-p2),2)
mtext(text=paste0('Q=',Q),side=3,cex=0.8)
col <- c(rep('white',n),rep('gray50',n))
legend('topleft','a)',bty='n',adj=c(1.5,0),cex=1.5)
plot(x,y,pch=21,axes=TRUE,ann=TRUE,col='gray50',
     bty='n',lwd=lwd,xlim=c(mx-dx/50,Mx+dx/50),
     ylim=c(my-dy/50,My+dy/50),bg=col,
     xaxs="i",yaxs="i",xlab='x',ylab='y')
plotframe(x,y)
s <- 0
lines(x=range(x),y=rep(s,2))
p1 <- sum(y[1:n] < s)/n
p2 <- sum(y[(n+1):(2*n)] < s)/n
Q <- signif(p1*(1-p1) + p2*(1-p2),2)
mtext(text=paste0('Q=',Q),side=3,cex=0.8)
col <- c(rep('white',n),rep('gray50',n))
legend('topleft','b)',bty='n',adj=c(1.5,0),cex=1.5)
plot(x,y,pch=21,axes=TRUE,ann=TRUE,col='gray50',
     bty='n',lwd=lwd,xlim=c(mx-dx/50,Mx+dx/50),
     ylim=c(my-dy/50,My+dy/50),bg=col,
     xaxs="i",yaxs="i",xlab='x',ylab='y')
plotframe(x,y)
s <- -4.902
lines(x=rep(s,2),y=range(y))
p1 <- sum(x[1:n] < s)/n
p2 <- sum(x[(n+1):(2*n)] < s)/n
Q <- signif(p1*(1-p1) + p2*(1-p2),2)
mtext(text=paste0('Q=',Q),side=3,cex=0.8)
legend('topleft','c)',bty='n',adj=c(1.5,0),cex=1.5)
dev.off()

cairo(file='../../figures/Q2.pdf',width=7,height=3.5)
pars(mfrow=c(1,2))
mx <- min(x)
Mx <- max(x)
my <- min(y)
My <- max(y)
plot(x,y,pch=21,axes=TRUE,ann=TRUE,col='gray50',
     bty='n',lwd=lwd,xlim=c(mx-dx/50,Mx+dx/50),
     ylim=c(my-dy/50,My+dy/50),bg=col,
     xaxs="i",yaxs="i",xlab='x',ylab='y',
     mar=c(2.5,2.3,0.5,0))
lines(x=range(x),y=rep(my,2))
lines(x=range(x),y=rep(My,2))
lines(x=rep(mx,2),y=range(y))
lines(x=rep(Mx,2),y=range(y))
lines(x=rep(-4.902,2),y=range(y))
lines(x=c(-4.903,Mx),y=rep(1.923,2))
my.control <- rpart.control(cp=0,minsplit=1,maxdepth=3)
tree <- rpart(groups ~ x + y, method="class", control=my.control)
plot(c(-1.6,1.1),c(-0.02,0.5),type='n',axes=FALSE,bty='n',ann=FALSE,
     mar=c(0,0,0.5,0))
lines(c(0,0),c(0.45333,0.5))
lines(c(-1,1),c(0.45333,0.45333))
lines(c(-1,-1),c(0.45333,0.12000))
lines(c(1,1),c(0.45333,0.12000))
lines(c(-1.5,-1.5),c(0.12,0))
lines(c(-0.5,-0.5),c(0.12,0))
lines(c(-1.5,-0.5),c(0.12,0.12))
text(x=0,y=0.43333,labels=expression('x'>='-4.902'),cex=0.8)
text(x=-1,y=0.10,labels=expression('y'>='1.923'),cex=0.8)
text(x=-1.5,y=-0.02,labels='54/2',cex=0.8)
text(x=-0.5,y=-0.02,labels='96/80',cex=0.8)
text(x=1.0,y=0.10,labels='0/68',cex=0.8)
dev.off()

cairo(file='../../figures/overfittedCARTscatter.pdf',width=7,height=3.5)
pars(mfrow=c(1,2),mar=c(2.5,2.3,0,0))
mx <- min(x)
Mx <- max(x)
my <- min(y)
My <- max(y)
plot(x,y,pch=21,axes=TRUE,ann=TRUE,col='gray50',
     bty='n',lwd=lwd,xlim=c(mx-dx/50,Mx+dx/50),
     ylim=c(my-dy/50,My+dy/50),bg=col,
     xaxs="i",yaxs="i",xlab='x',ylab='y')
my.control <- rpart.control(xval=10, cp=0, minsplit=1)
set.seed(5)
tree.unpruned <- rpart(groups ~ x + y, method="class", control=my.control)
lines(x=range(x),y=rep(my,2))
lines(x=range(x),y=rep(My,2))
lines(x=rep(mx,2),y=range(y))
lines(x=rep(Mx,2),y=range(y))
lines(x=rep(-4.902,2),y=range(y))
lines(x=c(-4.903,Mx),y=rep(1.923,2))
lines(x=rep(-3.376,2),y=c(1.923,My))
lines(x=c(-3.376,Mx),y=rep(2.236,2))
lines(x=c(-3.376,Mx),y=rep(2.175,2))
lines(x=c(-4.902,Mx),y=rep(-1.639,2))
lines(x=rep(-3.632,2),y=c(my,-1.639))
lines(x=rep(1.1,2),y=c(my,-1.639))
lines(x=rep(1.157,2),y=c(my,-1.639))
lines(x=rep(1.35,2),y=c(my,-1.639))
lines(x=rep(1.265,2),y=c(my,-1.639))
lines(x=c(1.35,Mx),y=rep(-2.314,2))
lines(x=c(1.35,Mx),y=rep(-2.229,2))
lines(x=rep(-0.5591,2),y=c(-1.639,1.923))
lines(x=rep(-4.081,2),y=c(-1.639,1.923))
lines(x=rep(-1.312,2),y=c(-1.639,1.923))
lines(x=c(-1.312,-0.5591),y=rep(-0.4715,2))
lines(x=c(-4.902,-4.081),y=rep(-0.0939,2))
lines(x=c(-4.902,-4.081),y=rep(-1.238,2))
lines(x=c(-4.902,-4.081),y=rep(1.775,2))
lines(x=rep(-4.613,2),y=c(-0.0939,1.775))
lines(x=rep(-4.332,2),y=c(-0.0939,1.775))
lines(x=c(-0.5591,Mx),y=rep(-1.409,2))
lines(x=c(-0.5591,Mx),y=rep(-1.523,2))
lines(x=c(-0.5591,Mx),y=rep(1.7,2))
lines(x=rep(1.04,2),y=c(1.7,1.923))
lines(x=rep(1.621,2),y=c(-1.409,1.7))
lines(x=rep(1.519,2),y=c(-1.409,1.7))
lines(x=c(-0.5591,1.519),y=rep(-1.077,2))
lines(x=rep(0.2677,2),y=c(-1.409,-1.077))
lines(x=rep(1.289,2),y=c(-1.077,1.7))
lines(x=rep(1.321,2),y=c(-1.077,1.7))
plot(tree.unpruned)
#text(tree.unpruned,use.n=F,xpd=NA)
dev.off()

cairo(file='../../figures/cvCART.pdf',width=3,height=3)
pars(mar=c(2.5,2.3,2.3,0.2),mgp=c(1.5,0.5,0))
plotcp(tree.unpruned,xpd=NA)
mtext('size of tree',side=3,line=1.5)
dev.off()

cairo(file='../../figures/optimalCART.pdf',width=7,height=3.5)
set.seed(5)
pars(mfrow=c(1,2),mar=c(2.5,2.3,0.5,0))
mx <- min(x)
Mx <- max(x)
my <- min(y)
My <- max(y)
plot(x,y,pch=21,axes=TRUE,ann=TRUE,col='gray50',
     bty='n',lwd=lwd,xlim=c(mx-dx/50,Mx+dx/50),
     ylim=c(my-dy/50,My+dy/50),bg=col,
     xaxs="i",yaxs="i",xlab='x',ylab='y')
lines(x=range(x),y=rep(my,2))
lines(x=range(x),y=rep(My,2))
lines(x=rep(mx,2),y=range(y))
lines(x=rep(Mx,2),y=range(y))
lines(x=rep(-4.902,2),y=range(y))
lines(x=c(-4.902,Mx),y=rep(1.923,2))
lines(x=c(-4.902,Mx),y=rep(-1.639,2))
lines(x=rep(-0.5591,2),y=c(-1.639,1.923))
tree <- rpart(groups ~ x + y, method="class")
plot(tree,xpd=NA,margin=0.05)
text(tree,use.n=T,xpd=NA,pos=1,offset=-0.19)
dev.off()

cairo(file='../../figures/irisCART.pdf',width=2.7,height=2.7)
pars(mar=c(1.5,1.1,1.1,1.1))
tree <- rpart(Species ~ ., data=IRIS, method="class")
plot(tree)
text(tree,use.n=TRUE,xpd=NA,pos=1,offset=-0.19)
dev.off()

if (FALSE){
    right <- (iris[,'Petal.Length']>=2.45)
    p <- par(mfrow=c(1,3))
    plot(iris[right,'Petal.Length'],iris[right,'Sepal.Length'])
    plot(iris[right,'Petal.Length'],iris[right,'Sepal.Width'])
    plot(iris[right,'Petal.Length'],iris[right,'Petal.Width'])
    par(p)
}

set.seed(1)
A <- runif(10)
B <- runif(10)
AB <- A/B
BA <- B/A
mAB <- mean(AB)
mBA <- mean(BA)
mBAinv <- 1/mBA
lAB <- log(A/B)
lBA <- log(B/A)
mlAB <- mean(lAB)
mlBA <- mean(lBA)
expmlAB <- exp(mlAB)
expmlBA <- exp(mlBA)
tab <- rbind(c(A,NA,NA),c(B,NA,NA),c(AB,mAB,NA),
             c(BA,mBA,mBAinv),c(lAB,mlAB,expmlAB),
             c(lBA,mlBA,expmlBA))
#signif(tab,2)

clr <- function(dat,inverse=FALSE){
    if (class(dat)%in%c('matrix','data.frame')){
        d <- dat
    } else {
        d <- matrix(dat,nrow=1)
    }
    if (inverse){
        edat <- exp(d)
        sdat <- rowSums(edat)
        out <- apply(edat,MARGIN=2,FUN='/',sdat)
    } else {
        lg <- rowMeans(log(d))
        out <- apply(log(d),MARGIN=2,FUN='-',lg)
    }
    as.matrix(out)
}
alr <- function(dat,inverse=FALSE){
    if (class(dat)%in%c('matrix','data.frame')){
        d <- dat
    } else {
        d <- matrix(dat,nrow=1)
    }
    if (inverse){
        num <- cbind(1,exp(d))
        den <- 1+rowSums(exp(d),na.rm=TRUE)
        out <- num/den        
    } else {
        out <- log(d[,-1])-log(d[,1])
    }
    as.matrix(out)
}
xyz2xy <- function(xyz){
    if (class(xyz)%in%c('matrix','data.frame')){
        n <- nrow(xyz)
        x <- xyz[,1]
        y <- xyz[,2]
        z <- xyz[,3]
    } else {
        n <- 1
        x <- xyz[1]
        y <- xyz[2]
        z <- xyz[3]
    }
    xy <- matrix(0,nrow=n,ncol=2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- sin(pi/3)*x/(x+y+z)
    return(xy)
}
ternary <- function(xyz=NULL,f=rep(1,3),labels=c('X','Y','Z'),
                    add=FALSE,type='p',...){
    if (class(xyz)%in%c('matrix','data.frame')){
        xyz <- as.matrix(xyz,nrow=1)
    }
    if (!add){
        corners <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,0,0))
        xy <- xyz2xy(corners)
        graphics::plot(xy,type='l',asp=1,axes=FALSE,
                       ann=FALSE,bty='n')
        position <- c(3,1,1)
        for (i in 1:3){
            if (f[i]==1) lab <- labels[i]
            else lab <- paste0(f[i],'x',labels[i])
            graphics::text(xy[i,,drop=FALSE],labels=lab,pos=position[i],xpd=NA)
        }
    }
    XYZ <- t(apply(matrix(xyz,ncol=3),MARGIN=1,FUN='*',f))
    XYZ <- apply(XYZ,MARGIN=2,FUN='/',rowSums(XYZ))
    if (!is.null(xyz)){
        if (type=='p'){
            graphics::points(xyz2xy(XYZ),...)
        } else if (type=='l'){
            graphics::lines(xyz2xy(XYZ),...)
        } else if (type=='t'){
            graphics::text(xyz2xy(XYZ),labels=rownames(xyz),...)
        }
    }
}

cairo(file='../../figures/ACNK.pdf',width=3.5,height=3)
pars(mar=c(1,1.5,1,1))
ACNK <- read.csv('ACNK.csv',header=TRUE,row.names=1,check.names=FALSE)
ternary(ACNK,type='p',
        labels=c(expression('Al'[2]*'O'[3]),
                 expression('CaO+Na'[2]*'O'),
                 expression('K'[2]*'O')))
dev.off()

cairo(file='../../figures/ACNKarithmeticmean.pdf',width=3.5,height=3)
pars(mar=c(1,1.5,1,1))
ACNK <- read.csv('ACNK.csv',header=TRUE,row.names=1,check.names=FALSE)
ternary(ACNK,type='p',
        labels=c(expression('Al'[2]*'O'[3]),
                 expression('CaO+Na'[2]*'O'),
                 expression('K'[2]*'O')),col='gray50')
mu <- colMeans(ACNK)
ternary(mu,add=TRUE,pch=22,bg='black')
dev.off()

cairo(file='../../figures/ACNKnaive.pdf',width=3.5,height=3)
pars(mar=c(1,1.5,1,1))
ACNK <- read.csv('ACNK.csv',header=TRUE,row.names=1,check.names=FALSE)
ternary(ACNK,type='p',
        labels=c(expression('Al'[2]*'O'[3]),
                 expression('CaO+Na'[2]*'O'),
                 expression('K'[2]*'O')),col='gray50')
sig <- apply(ACNK,MARGIN=2,FUN='sd')
LL <- mu - 2*sig
UL <- mu + 2*sig
ternary.polygon(LL,UL,col='black')
dev.off()

cairo(file='../../figures/alr.pdf',width=10,height=3.5)
m <- rbind(c(0.05,0.4,0.0,1.00),
           c(0.70,1.00,0.15,0.95),
           c(0.00,1.00,0.00,1.00))
split.screen(m)
ACNK <- read.csv('ACNK.csv',header=TRUE,row.names=1,check.names=FALSE)
uv <- alr(ACNK)
mu <- colMeans(uv)
s <- cov(uv)
ell <- IsoplotR::ellipse(x=mu[1],y=mu[2],covmat=s,alpha=0.05,n=100)
gmu <- alr(matrix(mu,ncol=2),inverse=TRUE)
gell <- alr(ell,inverse=TRUE)
screen(1); par(mar=rep(0,4))
ternary(ACNK,type='p',
        labels=c(expression('z=Al'[2]*'O'[3]),
                 expression('x=CaO+Na'[2]*'O'),
                 expression('y=K'[2]*'O')),col='gray50')
ternary(gmu,add=TRUE,type='p',pch=22,bg='black')
ternary(gell,add=TRUE,type='l')
screen(2,new=FALSE); par(mar=rep(0,4))
plot(ell,type='n',bty='n',axes=FALSE)
axis(side=1,mgp=c(1.5,0.5,0))
mtext(expression('u=ln[(CaO+Na'[2]*'O)/Al'[2]*'O'[3]*']'),side=1,line=1.5)
axis(side=2,mgp=c(1.5,0.5,0))
mtext(expression('v=ln[K'[2]*'O/Al'[2]*'O'[3]*']'),side=2,line=1.2)
points(uv,col='gray50')
polygon(ell)
points(mu[1],mu[2],pch=22,bg='black')
screen(3,new=FALSE); par(mar=rep(0,4))
plot(c(0,1),c(0,1),type='n',bty='n',ann=FALSE,axes=FALSE)
arrows(0.4,0.65,0.58,0.65,length=0.1)
text(0.49,0.65,labels='logratio transformation',pos=1,offset=-1)
arrows(0.58,0.6,0.4,0.6,length=0.1)
text(0.49,0.6,labels='inverse logratio transformation',pos=1,offset=1)
close.screen(3,all.screens=TRUE)
dev.off()

cairo(file='../../figures/abc.pdf',width=2.8,height=2.5)
pars(mar=rep(1,4))
xyz <- rbind(c(0.03,99.88,0.09),
             c(70.54,25.95,3.51),
             c(72.14,26.54,1.32))
colnames(xyz) <- c('a','b','c')
rownames(xyz) <- 1:3
ternary(xyz,type='p',labels=colnames(xyz),cex=2,pch=21,bg='white',f=c(1,1,5))
ternary(xyz,type='t',add=TRUE,f=c(1,1,5))
dev.off()

cairo(file='../../figures/alrPCA.pdf',width=3,height=3)
pars(mar=c(2.5,2.3,1.8,1.8))
Xa <- cbind(log(xyz[,'a']/xyz[,'c']),log(xyz[,'b']/xyz[,'c']))
colnames(Xa) <- c('     ln(a/c)','ln(b/c)     ')
pc <- prcomp(Xa)
biplot(pc,asp=1.4,col=c('gray50','black'))
dev.off()

cairo(file='../../figures/clrPCA.pdf',width=3,height=3)
pars(mar=c(2.5,2.3,1.8,1.8))
g <- exp(rowMeans(log(xyz)))
Xc <- cbind(log(xyz[,'a']/g),log(xyz[,'b']/g),log(xyz[,'c']/g))
colnames(Xc) <- c('','','')
pc <- prcomp(Xc)
biplot(pc,arrow.len=0.08,asp=0.9,col=c('gray50','black'))
text(x=0,y=-1,labels='c',col='black',cex=1.1)
text(x=-4,y=0,labels='b',col='black',cex=1.1)
text(x=4,y=0,labels='a',col='black',cex=1.1)
dev.off()

cairo(file='../../figures/majorPCA.pdf',width=4,height=4)
pars(mar=c(2.5,2.3,1.8,1.8))
Major <- read.csv(file="Major.csv",header=TRUE,row.names=1)
cMajor <- log(Major) - rowMeans(log(Major)) %*% matrix(1,1,ncol(Major))
pc <- prcomp(cMajor)
biplot(pc,asp=1.1,col=c('gray50','black'))
dev.off()

cairo(file='../../figures/AFM.pdf',width=3,height=3)
pars(mar=c(0,0,0,0))
afm <- read.csv('AFM.csv',header=TRUE,check.names=FALSE)
AFFINITY <- afm[,1]
th <- (AFFINITY=='th')
ca <- (AFFINITY=='ca')
col <- rep('black',length(AFFINITY))
col[th] <- 'white'
ternary(afm[,-1],type='p',pch=21,bg=col,labels=c('F','A','M'),cex=0.6)
dev.off()

cairo(file='../../figures/lrAFM.pdf',width=3,height=3)
pars()
lr <- alr(afm[,-1])
plot(lr,pch=21,bg=col,xlab='ln(A/F)',ylab='ln(M/F)',cex=0.6)
dev.off()

cairo(file='../../figures/LDAAFM.pdf',width=6,height=3)
pars(mfrow=c(1,2))
col <- rep('gray50',length(AFFINITY))
col[th] <- 'white'
plot(lr,pch=21,cex=0.6,bg=col,col='gray50',xlab='ln(A/F)',ylab='ln(A/F)')
ldafm <- construct_DA(lr,AFFINITY=afm[,'affinity'],padding=2,
                      quadratic=FALSE,plot=FALSE)
lines(ldafm$contours[[1]])
legend('topleft',legend='a)',bty='n',adj=c(2,0))
p <- par(mar=rep(0,4))
ternary(afm[,-1],pch=21,bg=col,cex=0.6,col='gray50',labels=c('F','A','M'))
terncont <- alr(ldafm$contours[[1]],inverse=TRUE)
ternary(terncont,type='l',add=TRUE)
legend('topleft',legend='b)',bty='n')
par(p)
dev.off()

cairo(file='../../figures/cath.pdf',width=6,height=3)
pars(mfrow=c(1,2))
p <- list()
p[[1]] <- list(uv0=c(-2,1),lA=1,lF=2,lM=3,lab='i',lty=1,srt=-36.5,pos=3,o=0.1)
p[[2]] <- list(uv0=c(-2,1),lA=1,lF=2,lM=5,lab='i',lty=2,srt=-36.5,pos=1,o=0.3)
p[[3]] <- list(uv0=c(1,2),lA=1,lF=2,lM=2,lab='ii',lty=3,srt=0,pos=3,o=0)
tt <- seq(from=0,to=5,length.out=100)
plot(c(-3,7),c(-15,5),type='n',xlab='ln(F/A)',ylab='ln(M/A)')
for (i in 1:length(p)){
    uv0 <- p[[i]]$uv0
    lF <- p[[i]]$lF
    lA <- p[[i]]$lA
    lM <- p[[i]]$lM
    u <- uv0[1] + (lF-lA)*tt
    v <- uv0[2] + (lF-lM)*tt
    lines(u,v,lty=p[[i]]$lty)
    points(uv0[1],uv0[2],pch=21,bg='white',cex=2)
    text(uv0[1],uv0[2],labels=p[[i]]$lab)
#    text(u[60],v[60],
#         labels=bquote(lambda[A]*'='*.(lA)*', '*
#                       lambda[F]*'='*.(lF)*', '*
#                       lambda[M]*'='*.(lM)),
#         pos=p[[i]]$pos,offset=p[[i]]$o,srt=p[[i]]$srt,col='gray50')
}
P <- par(mar=c(0.5,0,0.5,0))
ternary(alr(lr,inverse=TRUE),type='n',labels=c('A','F','M'))
for (i in 1:length(p)){
    uv0 <- p[[i]]$uv0
    lF <- p[[i]]$lF
    lA <- p[[i]]$lA
    lM <- p[[i]]$lM
    u <- uv0[1] + (lF-lA)*tt
    v <- uv0[2] + (lF-lM)*tt
    ternary(alr(cbind(u,v),inverse=TRUE),type='l',add=TRUE,lty=p[[i]]$lty)
    AFM0 <- alr(matrix(uv0,nrow=1),inverse=TRUE)
    rownames(AFM0) <- p[[i]]$lab
    ternary(AFM0,type='p',add=TRUE,pch=21,bg='white',cex=2)
    ternary(AFM0,type='t',add=TRUE)
}
par(P)
dev.off()

cairo(file='../../figures/PCAiris.pdf',width=3,height=3)
p <- par(mar=c(2.5,2.5,1.5,1.5),mgp=c(1.5,0.5,0))
pc <- prcomp(iris[,-5])
biplot(pc,col=c('gray60','black'),expand=2,
       xlim=c(-0.15,0.45),ylim=c(-0.25,0.23))
par(p)
dev.off()

cairo(file='../../figures/elbow.pdf',width=3,height=3)
pars()
K <- 10            # maximum number of clusters to evaluate
ss <- rep(0,K)     # initialise the vector with the sums of squares
for (k in 1:K){    # loop through all the k values
  fit <- kmeans(iris[,-5],centers=k) # fit the k means
  ss[k] <- fit$tot.withinss          # extract the sum of squares
}
plot(x=1:K,y=ss,type='b')            # plot as both lines and points
dev.off()

cairo(file='../../figures/DZtree.pdf',width=3,height=2)
pars(mar=c(0,1.5,1,0))
np <- length(DZ)                    # there are 5 samples in DZ
KS <- matrix(0,np,np)               # initialise the matrix of K-S values
snames <- names(DZ)                 # get the sample names
rownames(KS) <- snames              # label the rows of the KS matrix
colnames(KS) <- snames              # label the columns of the KS matrix
for (i in 1:np){                    # loop through the rows
  for (j in 1:np){                  # loop through the columns
    if (i!=j) KS[i,j] <- ks.test(DZ[[i]],DZ[[j]])$statistic
  }
}
d <- as.dist(KS)
tree <- hclust(d)
plot(tree,main='',sub='',xlab='',ylab='')
dev.off()

cairo(file='../../figures/testPCA.pdf',width=4,height=4)
pars(mar=c(2.5,2.5,1.5,1.5))
data(test,package='geostats')
lrdat <- clr(test[,-1])
pc <- prcomp(lrdat)
biplot(pc,col=c('grey60','black'),xlabs=test[,1],
       xlim=c(-0.31,0.2),ylim=c(-0.15,0.37))
dev.off()

pebbles <- c(44,51,79,65,27,31,4,355,22,352,287,
             7,287,339,0,276,342,355,334,296,7,
             17,351,349,37,339,40,324,325,334)

plot.circ <- function(angles,degrees=FALSE,tl=0.1,...){
    plot(x=c(-1,1),y=c(-1,1),type='n',axes=FALSE,
         ann=FALSE,asp=1,bty='n')
    symbols(0,0,circles=1,add=TRUE,inches=FALSE)
    ticks.circ(angles,degrees=degrees,tl=tl,...)
    lines(c(0,0),c(1,1+tl))
    lines(c(0,0),-c(1,1+tl))
    lines(c(1,1+tl),c(0,0))
    lines(-c(1,1+tl),c(0,0))
    text(0,1+tl,labels='0',pos=3,xpd=NA,offset=0.1)
    text(1+tl,0,labels='90',pos=4,xpd=NA,offset=0.1)
    text(0,-1-tl,labels='180',pos=1,xpd=NA,offset=0.15)
    text(-1-tl,0,labels='270',pos=2,xpd=NA,offset=0.1)
}

ticks.circ <- function(angles,degrees=FALSE,tl=0.1,...){
    if (degrees) rads <- angles*pi/180
    else rads <- angles
    x1 <- sin(rads)
    y1 <- cos(rads)
    x2 <- x1*(1-tl)
    y2 <- y1*(1-tl)
    matlines(x=rbind(x1,x2),y=rbind(y1,y2),
             lty=1,col='black',...)
}

points.circ <- function(angles,degrees=FALSE,tl=0.1,...){
    if (degrees) rads <- angles*pi/180
    else rads <- angles
    points(x=sin(rads),y=cos(rads),...)
}

addarrow <- function(xy0=c(0,0),a=0,plot=TRUE){
    x1 <- xy0[1] + cos(a*pi/180)
    y1 <- xy0[2] + sin(a*pi/180)
    if (plot)
        arrows(x0=xy0[1],y0=xy0[2],x1=x1,y1=y1,length=0.1)
    invisible(c(x1,y1))
}

cairo(file='../../figures/circle1.pdf',width=2,height=2)
pars(mar=c(0,1.5,0,1.1))
plot.circ(pebbles,degrees=TRUE)
dev.off()

cairo(file='../../figures/circle2.pdf',width=2,height=2)
pars(mar=c(0,1.5,0,1.1))
plot.circ(pebbles,degrees=TRUE)
addarrow(a=90-mean(pebbles))
dev.off()

cairo(file='../../figures/circle3.pdf',width=2,height=2)
pars(mar=c(0,1.5,0,1.1))
rad <- pebbles*pi/180 # convert to radians
ss <- sum(sin(rad))
# sum of the sines
sc <- sum(cos(rad))
# sum of the cosines
md <- atan(ss/sc)
# mean direction
deg <- md*180/pi
# convert radians to degrees
plot.circ(pebbles,degrees=TRUE)
addarrow(a=90-deg)
dev.off()

cairo(file='../../figures/vectorsum.pdf',width=5,height=2.5)
pars(mar=rep(0,4))
plot(x=c(-1,3.5),y=c(-1,1),type='n',axes=FALSE,ann=FALSE,bty='n',asp=1)
symbols(0,0,circles=1,add=TRUE,inches=FALSE,xpd=NA)
xy1 <- addarrow(a=45)
xy2 <- addarrow(a=30)
xy3 <- addarrow(a=0)
xy4 <- addarrow(a=-30)
xy5 <- addarrow(xy0=xy4,a=0)
xy6 <- addarrow(xy0=xy5,a=30)
xy7 <- addarrow(xy0=xy6,a=45)
lines(x=c(xy3[1],xy5[1]),y=c(xy3[2],xy5[2]),lty=3)
lines(x=c(0,xy5[1]),y=c(0,xy5[2]),lty=3)
lines(x=c(xy2[1],xy6[1]),y=c(xy2[2],xy6[2]),lty=3)
lines(x=c(0,xy6[1]),y=c(0,xy6[2]),lty=3)
lines(x=c(xy1[1],xy7[1]),y=c(xy1[2],xy7[2]),lty=3)
arrows(x0=0,x1=xy7[1],y0=0,y1=xy7[2],lty=2,lwd=1,length=0)
arrows(x0=0,x1=xy7[1]/4,y0=0,y1=xy7[2]/4,
       lty=1,lwd=2,col='black',length=0.1)
text(xy7[1]/4,xy7[2]/4,labels=expression(bar(R)),pos=3,offset=0.2)
a <- seq(from=0,to=atan(xy7[2]/xy7[1]),length.out=10)
x <- 1.05*cos(a)
y <- 1.05*sin(a)
lines(x,y)
text(x[4],y[4],labels=expression(bar(theta)),font=5,pos=4,offset=0.1)
dev.off()

cairo(file='../../figures/lowconcentration.pdf',width=2,height=2)
pars(mar=rep(0,4))
plot(x=c(-1,1),y=c(-1,1),type='n',axes=FALSE,ann=FALSE,bty='n',asp=1)
symbols(0,0,circles=1,add=TRUE,inches=FALSE,xpd=NA)
xy1 <- addarrow(a=230)
xy2 <- addarrow(a=30)
xy3 <- addarrow(a=0)
xy4 <- addarrow(a=-170)
xy5 <- addarrow(xy0=xy4,a=0,plot=FALSE)
xy6 <- addarrow(xy0=xy5,a=30,plot=FALSE)
xy7 <- addarrow(xy0=xy6,a=230,plot=FALSE)
arrows(x0=0,x1=xy7[1]/4,y0=0,y1=xy7[2]/4,lty=1,lwd=2,col='black',length=0.1)
text(xy7[1]/4,xy7[2]/4,labels=expression(bar(R)),pos=1,offset=0.3)
dev.off()

vonMises <- function(angle,mu=0,kappa=1){
    num <- exp(kappa*cos(angle-mu))
    den <- 2*pi*besselI(kappa,nu=0)
    num/den
}

cairo(file='../../figures/vonMises.pdf',width=8,height=2)
pars(mar=rep(0,4))
plot(x=c(-1,11),y=c(-1.1,1.1),type='n',axes=FALSE,ann=FALSE,bty='n',asp=1)
a <- seq(from=-pi,to=pi,length.out=200)
d <- vonMises(angle=a,mu=pi/2,kappa=0)
x0 <- 0
symbols(x=x0,y=0,circles=1,add=TRUE,inches=FALSE,xpd=NA,fg='grey50')
text(x=x0,y=0,labels=expression(kappa*'=0'))
lines(x=x0+(1+d)*cos(a),y=(1+d)*sin(a),xpd=NA)
d <- vonMises(angle=a,mu=-pi/2,kappa=1)
x0 <- 3
symbols(x=x0,y=0,circles=1,add=TRUE,inches=FALSE,xpd=NA,fg='grey50')
text(x=x0,y=0,labels=expression(mu*'='*pi*','~kappa*'=1'))
lines(x=x0+(1+d)*cos(a),y=(1+d)*sin(a),xpd=NA)
d <- vonMises(angle=a,mu=pi/4,kappa=2)
x0 <- 6
symbols(x=x0,y=0,circles=1,add=TRUE,inches=FALSE,xpd=NA,fg='grey50')
text(x=x0,y=0,labels=expression(mu*'='*pi*'/4,'~kappa*'=2'))
lines(x=x0+(1+d)*cos(a),y=(1+d)*sin(a),xpd=NA)
d <- vonMises(angle=a,mu=0,kappa=10)
x0 <- 9
symbols(x=x0,y=0,circles=1,add=TRUE,inches=FALSE,xpd=NA,fg='grey50')
text(x=x0,y=0,labels=expression(mu*'='*pi*'/2,'~kappa*'=10'))
lines(x=x0+(1+d)*cos(a),y=(1+d)*sin(a),xpd=NA)
dev.off()

cairo(file='../../figures/R2K.pdf',width=2.5,height=2.5)
pars()
R<-seq(from=0.01,to=0.99,by=0.01)
plot(x=R,y=R*(2-R^2)/(1-R^2),type='l',
     xlab=expression(bar(R)),ylab=expression(kappa))
dev.off()

cairo(file='../../figures/kappastriations.pdf',width=3.5,height=2.5)
mu <- geostats::meanangle(pebbles,degrees=TRUE)
Rp <- geostats::Rbar(pebbles,degrees=TRUE)
Kp <- geostats::Rbar2kappa(Rp)
a <- seq(from=-pi,to=pi,length.out=200)
f <- vonMises(angle=a,mu=mu*pi/180,kappa=Kp)
pars(mar=c(3,3,0.5,0.5))
plot(a,f,type='l',bty='n',xaxt='n')
axis(side=1,at=c(-pi,-pi/2,0,pi/2,pi),
     labels=c(expression(-pi),expression(-pi/2),0,expression(pi/2),expression(pi)))
A <- atan(sin(pebbles*pi/180)/cos(pebbles*pi/180))
rug(A,ticksize=0.1)
if (FALSE){
    plot(x=c(-1.1,1.1),y=c(-1.1,1.1),type='n',axes=FALSE,ann=FALSE,bty='n',asp=1)
    plot.circ(pebbles,degrees=TRUE)
    text(x=x0,y=0,labels=expression(kappa*'=0'))
    lines(x=x0+(1+f)*cos(a+pi/2),y=(1+f)*sin(a+pi/2),xpd=NA)
}
dev.off()

cairo(file='../../figures/wulffschmidt.pdf',width=7,height=3.5)
pars(mfrow=c(1,2),mar=rep(1,4))
geostats::stereonet(wulff=TRUE,show.grid=TRUE)
geostats::stereonet(trd=c(0,0,90,180,270),plg=c(90,10,10,10,10),
                    coneAngle=rep(10,5),option=4,
                    degrees=TRUE,add=TRUE,wulff=TRUE,lwd=1.5)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
geostats::stereonet(wulff=FALSE,show.grid=TRUE)
geostats::stereonet(trd=c(0,0,90,180,270),plg=c(90,10,10,10,10),
                    coneAngle=rep(10,5),option=4,
                    degrees=TRUE,add=TRUE,wulff=FALSE,lwd=1.5)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

cairo(file='../../figures/Africa.pdf',width=7,height=3.5)
pars(mfrow=c(1,2),mar=rep(1,4))
Africa <- read.csv('~/Documents/Programming/R/geostats/build/notes/Africa.csv',
                   header=TRUE)
geostats::stereonet(trd=Africa$lon,plg=Africa$lat,option=3,
                    degrees=TRUE,wulff=TRUE,type='l',lty=1.5)
legend('topleft',legend='a)',bty='n',cex=1.2,adj=c(2,0))
geostats::stereonet(trd=Africa$lon,plg=Africa$lat,option=3,
                    degrees=TRUE,wulff=FALSE,type='l',lty=1.5)
legend('topleft',legend='b)',bty='n',cex=1.2,adj=c(2,0))
dev.off()

xyz2SD <- function(xyz){
    D <- asin(xyz[,3])
    S <- sign(xyz[,2])*acos(xyz[,1]/cos(D))
    cbind(S,D)*180/pi
}
xyz2AD <- function(xyz){
    D <- asin(xyz[,3])
    A <- acos(xyz[,2]/cos(D))
    cbind(A,D)*180/pi
}
if (FALSE){
    library(Rfast)
    xyz <- c(10,10,5); mu <- xyz/sqrt(sum(xyz^2))
    palaeomag <- rvmf(n=10,mu=xyz,k=200)
    AD <- xyz2AD(palaeomag)
    A <- AD[,1]
    D1 <- AD[,2]
} else {
    A <- c(47.9,46.3,44.7,50.9,56.4,42.6,44.9,41.5,47.9,39.6)
    D1 <- c(28.6,20.1,15.6,18.1,17.5,28.7,12.2,24.5,20.6,15.0)
}
cairo(file='../../figures/palaeomag.pdf',width=2.5,height=2.5)
pars(mar=c(1.2,1,1,1))
geostats::stereonet(trd=A,plg=D1,option=1,degrees=TRUE,
                    show.grid=FALSE,wulff=FALSE,pch=19,cex=0.7)
dev.off()

if (FALSE){
    xyz <- c(-10,-10,10); mu <- xyz/sqrt(sum(xyz^2))
    fault <- rvmf(n=10,mu=xyz,k=100)
    SD <- xyz2SD(fault)
    S <- SD[,1]
    D2 <- SD[,2]    
} else {
    S <- c(226,220,223,222,233,227,234,229,227,224)
    D2 <- c(28.4,35.3,41.0,39.6,48.3,34.7,34.5,36.0,34.2,28.7)
}
cairo(file='../../figures/fault.pdf',width=2.5,height=2.5)
pars(mar=c(1.2,1,1,1))
geostats::stereonet(trd=S,plg=D2,option=2,degrees=TRUE,
                    show.grid=FALSE,wulff=TRUE,pch=21,cex=0.7,bg='white')
dev.off()

cairo(file='../../figures/sphericalmean.pdf',width=5,height=2.5)
pars(mfrow=c(1,2),mar=c(1.2,1,1,1))
meanpalaeomag <- geostats::meanangle(trd=A,plg=D1,option=1,degrees=TRUE)
geostats::stereonet(trd=A,plg=D1,option=1,degrees=TRUE,
                    show.grid=FALSE,wulff=FALSE,pch=19,cex=0.7,col='grey70')
geostats::stereonet(trd=meanpalaeomag[1],plg=meanpalaeomag[2],
                    option=1,degrees=TRUE,add=TRUE,wulff=FALSE,
                    pch=15,cex=0.7,col='black')
legend('topleft',legend='a)',adj=c(2,0),bty='n')
meanfault <- geostats::meanangle(trd=S,plg=D2,option=2,degrees=TRUE)
geostats::stereonet(trd=S,plg=D2,option=2,degrees=TRUE,
                    show.grid=FALSE,wulff=TRUE,pch=19,cex=0.7,col='grey70')
geostats::stereonet(trd=meanfault[1],plg=meanfault[2],
                    option=2,degrees=TRUE,add=TRUE,wulff=TRUE,
                    pch=15,cex=0.7,col='black',lwd=1.5)
legend('topleft',legend='b)',adj=c(2,0),bty='n')
dev.off()

data('meuse',package='geostats')
X <- meuse$x
Y <- meuse$y
Z <- log(meuse$zinc)

cairo(file='../../figures/meusepoints.pdf',width=5,height=4)
pars(mar=c(1.5,3,2,0))
geostats::colourplot(X=X,Y=Y,Z=Z,colspec=grey.colors,
                     key.title=title('ln[Zn]'),cex=0.7,
                     extra={text(179850,331650,labels='?')})
dev.off()

if (FALSE){
    N <- length(X)
    dXY <- as.matrix(dist(cbind(X,Y)))
    dZ <- as.matrix(dist(cbind(Z,Z)))
    d <- matrix(dXY[upper.tri(dXY,diag=TRUE)],N,N) +
        matrix(dZ[lower.tri(dZ,diag=TRUE)],N,N)
    signif(d[c(1:10,N),c(1:10,N)],2)
}

cairo(file='../../figures/semivariogram.pdf',width=3,height=3)
pars(mgp=c(1.2,0.5,0))
geostats::semivariogram(x=X,y=Y,z=Z,plot=TRUE,fit=FALSE)
dev.off()

cairo(file='../../figures/semivariogramfit.pdf',width=7,height=1.75)
pars(mfrow=c(1,4),mgp=c(1.2,0.5,0))
geostats::semivariogram(x=X,y=Y,z=Z,plot=TRUE,fit=TRUE,model='exponential')
legend('topleft',legend='a)',bty='n',adj=c(2,0),cex=1.2)
geostats::semivariogram(x=X,y=Y,z=Z,plot=TRUE,fit=TRUE,model='gaussian')
legend('topleft',legend='b)',bty='n',adj=c(2,0),cex=1.2)
geostats::semivariogram(x=X,y=Y,z=Z,plot=TRUE,fit=TRUE,model='linear')
legend('topleft',legend='c)',bty='n',adj=c(2,0),cex=1.2)
svm <- geostats::semivariogram(x=X,y=Y,z=Z,plot=TRUE,fit=TRUE,model='spherical')
legend('topleft',legend='d)',bty='n',adj=c(2,0),cex=1.2)
dev.off()

cairo(file='../../figures/snr.pdf',width=3,height=3)
pars(mar=c(2.5,2.3,0.8,0.25),mgp=c(1.2,0.5,0))
plot(x=c(0,500,1500),y=c(0.1,0.6,0.6),type='l',lwd=3,
     xlab='lag',ylab=expression(gamma(lag)),
     xlim=c(-50,1550),ylim=c(0,0.62),bty='n')
arrows(x0=-50,y0=0,x1=-50,y1=0.1,length=0.05,angle=90,code=3)
arrows(x0=1550,y0=0,x1=1550,y1=0.6,length=0.05,angle=90,code=3)
arrows(x0=0,y0=0.61,x1=500,y1=0.61,length=0.05,angle=90,code=3)
lines(x=c(-50,0),y=c(0.1,0.1),lty=3)
lines(x=c(-50,1550),y=c(0,0),lty=3)
text(x=-50,y=0.05,labels=expression(nugget~(c[n])),pos=4)
lines(x=c(0,0),y=c(0.1,0.61),lty=3)
lines(x=c(500,500),y=c(0.6,0.61),lty=3)
text(x=250,y=0.61,labels=expression(range~(c[r])),pos=3,xpd=NA)
lines(x=c(1500,1550),y=c(0.6,0.6),lty=3)
text(x=1550,y=0.3,labels=expression(sill~(c[s])),pos=2)
dev.off()

#png(file='../../figures/meusecontour.png',type='cairo',
#    family="serif",width=1000,height=800)
cairo(file='../../figures/meusecontour.pdf',width=5,height=4)
pars(mar=c(1.5,3,2,0))
dX <- (max(X)-min(X))/10
dY <- (max(Y)-min(Y))/10
xi <- seq(from=min(X)-dX,to=max(X)+dX,length.out=50)
yi <- seq(from=min(Y)-dY,to=max(Y)+dY,length.out=50)
zi <- geostats::kriging(x=X,y=Y,z=Z,svm=svm,xi=xi,yi=yi,grid=TRUE)
mZ <- min(zi,na.rm=TRUE)
MZ <- max(zi,na.rm=TRUE)
geostats::colourplot(x=xi,y=yi,z=zi,X=X,Y=Y,Z=Z,colspec=grey.colors,
                     key.title=title('ln[Zn]'),cex=0.7)
dev.off()

if (TRUE){
    znew <- geostats::kriging(x=X,y=Y,z=Z,svm=svm,
                              xi=179850,yi=331650,err=FALSE)
    sznew <- geostats::kriging(x=X,y=Y,z=Z,svm=svm,
                               xi=179850,yi=331650,err=TRUE)
}

cairo(file='../../figures/meusecontourerr.pdf',width=5,height=4)
pars(mar=c(1.5,3.5,2,0))
dX <- (max(X)-min(X))/10
dY <- (max(Y)-min(Y))/10
xi <- seq(from=min(X)-dX,to=max(X)+dX,length.out=50)
yi <- seq(from=min(Y)-dY,to=max(Y)+dY,length.out=50)
zi <- sqrt(geostats::kriging(x=X,y=Y,z=Z,svm=svm,xi=xi,yi=yi,grid=TRUE,err=TRUE))
mZ <- min(zi,na.rm=TRUE)
MZ <- max(zi,na.rm=TRUE)
geostats::colourplot(x=xi,y=yi,z=zi,colspec=grey.colors,cex=0.7,
                     key.title=title("s(Zn)/Zn",cex.main=1.04))
dev.off()

set.seed(8)
nn <- 150
mu1 <- c(10,10)
mu2 <- c(0,-10)
mu3 <- c(-10,10)
E <- diag(2)*30
xlim <- c(-25,25)
ylim <- c(-25,25)
X <- xlim[1] + runif(nn)*diff(xlim)
Y <- ylim[1] + runif(nn)*diff(ylim)
XY <- cbind(X,Y)
dolog <- FALSE
Z <- 400+1e5*(
    1.2*mvtnorm::dmvnorm(XY,mu1,E,log=dolog) +
    mvtnorm::dmvnorm(XY,mu2,E,log=dolog) +
    mvtnorm::dmvnorm(XY,mu3,E,log=dolog)*0.8
)
hills <- cbind(X,Y,Z)
colnames(hills) <- c('X','Y','Z')
cairo(file='../../figures/hills.pdf',width=8,height=2)
pars(mfrow=c(1,4),mar=c(2.5,2.3,0.7,0.25))
normodel <- geostats::semivariogram(x=X,y=Y,z=Z,model='gaussian',plot=TRUE)
legend('topleft',legend='a)',bty='n',adj=c(2,0),cex=1.2)
expmodel <- geostats::semivariogram(x=X,y=Y,z=Z,model='exponential',plot=TRUE)
legend('topleft',legend='b)',bty='n',adj=c(2,0),cex=1.2)
linmodel <- geostats::semivariogram(x=X,y=Y,z=Z,model='linear',plot=TRUE)
legend('topleft',legend='c)',bty='n',adj=c(2,0),cex=1.2)
sphmodel <- geostats::semivariogram(x=X,y=Y,z=Z,model='spherical',plot=TRUE)
legend('topleft',legend='d)',bty='n',adj=c(2,0),cex=1.2)
dev.off()

cairo(file='../../figures/sphvsexp.pdf',width=5,height=4)
pars(mar=c(1.5,3.5,2,0))
X <- meuse$x
Y <- meuse$y
Z <- log(meuse$zinc)
sphmodel <- geostats::semivariogram(x=X,y=Y,z=Z,plot=FALSE,model='spherical')
expmodel <- geostats::semivariogram(x=X,y=Y,z=Z,plot=FALSE,model='exponential')
xi <- seq(from=min(X),to=max(X),length.out=50)
yi <- seq(from=min(Y),to=max(Y),length.out=50)
zsph <- geostats::kriging(x=X,y=Y,z=Z,xi=xi,yi=yi,grid=TRUE,svm=sphmodel)
zexp <- geostats::kriging(x=X,y=Y,z=Z,xi=xi,yi=yi,grid=TRUE,svm=expmodel)
geostats::colourplot(x=xi,y=yi,z=exp(zsph-zexp),colspec=grey.colors,
                     extra={points(X,Y,pch=21,cex=0.7,bg='white')},
                     key.title=title(expression(Zn[sph]/Zn[exp]),cex.main=1.04))
dev.off()

cairo(file='../../slides/locationdispersionshape.pdf',width=4,height=2)
pars(mfrow=c(2,3),mgp=c(1.2,0.5,0))
x <- seq(from=-10,to=10,length.out=200)
y1 <- dnorm(x,mean=-2,sd=1)
y2 <- dnorm(x,mean=2,sd=1)
y3 <- dnorm(x,mean=0,sd=1)
y5 <- dlnorm(x,meanlog=1,sdlog=1/2)
y4 <- dnorm(x,mean=0,sd=2)
y6 <- dlnorm(x,meanlog=1,sdlog=1/2)
plot(x,y1,type='l',bty='n',ylim=c(0,0.5),ylab=expression('y'[1]))
plot(x,y3,type='l',bty='n',ylim=c(0,0.5),ylab=expression('y'[3]))
plot(x[100:200],y5[100:200],type='l',bty='n',ylim=c(0,0.5),
     ylab=expression('y'[5]),xlab='x')
plot(x,y2,type='l',bty='n',ylim=c(0,0.5),ylab=expression('y'[2]))
plot(x,y4,type='l',bty='n',ylim=c(0,0.5),ylab=expression('y'[4]))
plot(20-x[100:200],y6[100:200],type='l',bty='n',ylim=c(0,0.5),
     ylab=expression('y'[6]),xlab='x')
dev.off()

cairo(file='../../slides/birthdays.pdf',width=6,height=4)
pars()
k <- 1:100
lp <- lfactorial(365) - lfactorial(365-k) - k*log(365)
P <- 1 - exp(lp)
plot(k,P,type='l')
dev.off()

cairo(file='../../slides/uniform.pdf',width=6,height=2)
pars(mfrow=c(1,2),mar=c(2.5,2.3,1.0,0.2))
plot(c(-2,-1,-1,1,1,2),y=c(0,0,1,1,0,0),type='l',xlab='x',ylab='f',axes=FALSE)
if (TRUE){
    axis(side=1,at=c(-2,-1,1,2),labels=c('','a','b',''))
} else {
    axis(side=1,at=c(-2,-1,0.1,0.5,1,2),labels=c('','a','c','d','b',''))
    polygon(c(0.1,0.1,0.5,0.5,0.1),c(0,1,1,0,0),col='black')
}
axis(side=2,at=c(0,1),labels=c(0,'1/(b-a)'))
plot(c(-2,-1,1,2),y=c(0,0,1,1),type='l',xlab='x',ylab='F',bty='n',axes=FALSE)
axis(side=1,at=c(-2,-1,1,2),labels=c('','a','b',''))
axis(side=2)
dev.off()

cairo(file='../../slides/averageheight.pdf',width=5,height=3)
pars(mgp=c(1.2,0.5,0))
x <- 10^seq(from=0,to=9.69897,length.out=100)
y <- 10/sqrt(x)
plot(x,y,type='l',xlab='N',ylab=expression("s["*bar(h)*"]"),log='x',bty='n')
dev.off()

cairo(file='../../slides/RbSr.pdf',width=4,height=4)
pars(mgp=c(1.25,0.5,0))
RS <- t(RbSr[c(1,3,2,4,5),])
RS[,1] <- RS[,1]/2
RS[,5] <- 1.5*RS[,5]
scatterplot(RS,fit=york(RS),fill=NA,
            xlab=expression(''^87*'Rb/'^86*'Sr'),
            ylab=expression(''^87*'Sr/'^86*'Sr'))
dev.off()

cairo(file='../../slides/MDS2D.pdf',width=2.5,height=2.5)
pars(mar=c(2.5,2.5,0.5,0.5))
m <- cbind(c(4.24,-2.12,-2.12),c(0,-0.71,0.71))
plot(m,type='n',xlab='dim 1',ylab='dim 2')
text(m,labels=1:3)
dev.off()

cairo(file='../../slides/ChinaKDEs.pdf',width=6,height=4)
pars()
class(DZ) <- 'detritals'
IsoplotR::kde(DZ,kde.col=NA,hist.col=NA,show.hist=FALSE,
              rug=TRUE,samebandwidth=TRUE,normalise=FALSE)
dev.off()

cairo(file='../../slides/vectorsumnoR.pdf',width=4,height=3)
pars(mar=rep(0,4))
plot(x=c(-1,3.5),y=c(-1,1),type='n',axes=FALSE,ann=FALSE,bty='n',asp=1)
symbols(0,0,circles=1,add=TRUE,inches=FALSE,xpd=NA)
xy1 <- addarrow(a=45)
xy2 <- addarrow(a=30)
xy3 <- addarrow(a=0)
xy4 <- addarrow(a=-30)
xy5 <- addarrow(xy0=xy4,a=0)
xy6 <- addarrow(xy0=xy5,a=30)
xy7 <- addarrow(xy0=xy6,a=45)
lines(x=c(xy3[1],xy5[1]),y=c(xy3[2],xy5[2]),lty=3)
lines(x=c(0,xy5[1]),y=c(0,xy5[2]),lty=3)
lines(x=c(xy2[1],xy6[1]),y=c(xy2[2],xy6[2]),lty=3)
lines(x=c(0,xy6[1]),y=c(0,xy6[2]),lty=3)
lines(x=c(xy1[1],xy7[1]),y=c(xy1[2],xy7[2]),lty=3)
arrows(x0=0,x1=xy7[1],y0=0,y1=xy7[2],lty=1,lwd=2,length=0.1)
a <- seq(from=0,to=atan(xy7[2]/xy7[1]),length.out=10)
x <- 1.05*cos(a)
y <- 1.05*sin(a)
lines(x,y)
text(x[4],y[4],labels=expression(bar(theta)),font=5,pos=4,offset=0.1)
dev.off()

cairo(file='../../slides/Rbar.pdf',width=2.5,height=2.5)
pars(mar=rep(0,4))
plot(x=c(-1,1),y=c(-1,1),type='n',axes=FALSE,ann=FALSE,bty='n',asp=1)
symbols(0,0,circles=1,add=TRUE,inches=FALSE,xpd=NA)
xy1 <- addarrow(a=45)
xy2 <- addarrow(a=30)
xy3 <- addarrow(a=0)
xy4 <- addarrow(a=-30)
arrows(x0=0,x1=xy7[1]/4,y0=0,y1=xy7[2]/4,
       lty=1,lwd=2,col='black',length=0.1)
text(xy7[1]/4,xy7[2]/4,labels=expression(bar(R)),pos=3,offset=0.2)
dev.off()

# PDFs and CDFs for quizzes:

plotpdf <- function(x,d){
    plot(x,d,type='l',bty='n',ann=FALSE,axes=FALSE)
    Axis(side=1,at=c(0,0.5,1))
}
plotcdf <- function(x,d){
    dx <- diff(x)
    foo <- cumsum(d*c(dx[1],dx))
    cdf <- foo/tail(foo,n=1)
    plot(x,cdf,type='l',bty='n',ann=FALSE,axes=FALSE)
    Axis(side=1,at=c(0,0.5,1))
}
quiz1question <- function(qn=1){
    dx <- 1e-4
    x <- list()
    pd <- list()
    x[[1]] <- seq(from=0,to=1,length.out=100)
    pd[[1]] <- rep(1/100,100)
    x[[2]] <- c(seq(from=0,to=0.5-dx,length.out=50),
                0.5,
                seq(from=0.5+dx,to=1,length.out=50))
    foo <- c(rep(0,50),1,rep(0,50))
    pd[[2]] <- foo/sum(foo)
    set.seed(1)
    x[[3]] <- seq(from=0,to=1,length.out=100)
    pd[[3]] <- 0.3*dnorm(x[[3]],mean=0.3,sd=0.05) +
        0.7* dnorm(x[[3]],mean=0.7,sd=0.05)
    x[[4]] <- x[[3]]
    pd[[4]] <- 0.7*dnorm(x[[3]],mean=0.3,sd=0.05) +
        0.3* dnorm(x[[3]],mean=0.7,sd=0.05)
    x[[5]] <- seq(from=0,to=1,length.out=100)
    pd[[5]] <- x[[5]]/sum(x[[5]])
    pars(mfrow=c(6,2),mar=c(2,1,0.5,0.5))
    plotpdf(x[[qn]],pd[[qn]])
    legend('topleft',legend='1)',bty='n',adj=c(3,-1),xpd=NA,cex=1.2)
    plotcdf(x[[3]],pd[[3]])
    legend('topleft',legend='a)',bty='n',adj=c(3,-1),xpd=NA,cex=1.2)
    plot.new()
    plotcdf(x[[1]],pd[[1]])
    legend('topleft',legend='b)',bty='n',adj=c(3,-1),xpd=NA,cex=1.2)
    plot.new()
    plotcdf(x[[5]],pd[[5]])
    legend('topleft',legend='c)',bty='n',adj=c(3,-1),xpd=NA,cex=1.2)
    plot.new()
    plotcdf(x[[2]],pd[[2]])
    legend('topleft',legend='d)',bty='n',adj=c(3,-1),xpd=NA,cex=1.2)
    plot.new()
    plotcdf(x[[4]],pd[[4]])
    legend('topleft',legend='e)',bty='n',adj=c(3,-1),xpd=NA,cex=1.2)
    plot.new()
    plot(x=c(0,1),y=c(1,1),type='l',bty='n',ann=FALSE,axes=FALSE)
    Axis(side=1,at=c(0,0.5,1))
    legend('topleft',legend='f)',bty='n',adj=c(3,-1),xpd=NA,cex=1.2)
}

for (i in 1:5){
    fn <- paste0('~/Desktop/q1q',i,'.png')
    png(filename=fn,width=400,height=800,pointsize=14)
    quiz1question(qn=i)
    dev.off()
}

png(filename='~/Desktop/boxplot.png',width=400,height=200,pointsize=14)
pars(mar=c(2,0.5,0.5,0.5))
dat <- c(seq(from=0,to=50,length.out=25),
         seq(from=51,to=100,length.out=25),
         seq(from=101,to=200,length.out=25),
         seq(from=201,to=400,length.out=25)
         )
boxplot(dat,horizontal=TRUE,xaxt='n')
axis(1,at=c(0,50,100,200,400))
dev.off()

png(filename='~/Desktop/ecdf.png',width=400,height=400,pointsize=14)
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

png(filename='~/Desktop/ci.png',width=400,height=400,pointsize=14)
pars(mar=c(3,3,0,0))
xx <- c(1,2,3,4)
yy <- c(20,22,19,21)
dy <- c(4,1,3,2)
plot(xx,yy,type='n',ann=FALSE,bty='n',xaxt='n',
     xlim=c(0.5,4.5),ylim=range(c(yy-dy,yy+dy)))
arrows(xx,yy-dy,xx,yy+dy,code=3,angle=90)
axis(1,at=1:4,labels=c('a','b','c','d'))
dev.off()

png(filename='~/Desktop/QQ.png',width=800,height=800,pointsize=14)
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

png(filename='~/Desktop/pt.png',width=400,height=300,pointsize=14)
x <- seq(from=-5,to=5,length.out=100)
y <- pt(x,df=3)
plot(x,y,type='l',xlab='t',ylab='F')
lines(range(x),rep(0.025,2),lty=2)
lines(range(x),rep(0.975,2),lty=2)
lines(rep(qt(0.025,df=3),2),c(0,1),lty=3)
lines(rep(qt(0.975,df=3),2),c(0,1),lty=3)
dev.off()

png(filename='~/Desktop/majorqq.png',width=1200,height=400,pointsize=14)
pars(mfrow=c(2,5),mar=c(3,3,1,0.5))
data(major,package='geostats')
for (cname in colnames(major)){
    qqnorm(major[,cname],main=cname)
    qqline(major[,cname])
}
dev.off()
