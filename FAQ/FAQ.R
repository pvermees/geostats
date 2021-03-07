graphics.off()

setwd('~/Documents/Programming/R/geostats/build/notes/')

pars <- function(mar=c(2.5,2.3,0.5,0.25),mgp=c(1.5,0.5,0),mfrow=c(1,1)){
    par(list(mar=mar,mgp=mgp,mfrow=mfrow))
}

cairo <- function(file,width,height,family="serif",pointsize=13,...){
    cairo_pdf(file=file,width=width,height=height,
              family=family,pointsize=pointsize,...)
}

cairo(file='dof.pdf',width=4,height=4)
set.seed(1)
pars()
ns <- 1000
nv <- 2
obs <- matrix(rnorm(nv*ns),nrow=ns,ncol=nv)
tstat <- rep(NA,ns)
for (i in 1:ns){
    tstat[i] <- sqrt(nv)*mean(obs[i,])/sd(obs[i,])
}
pred1 <- qt(seq(from=0,to=1,length.out=ns),df=nv-1)
pred2 <- qt(seq(from=0,to=1,length.out=ns),df=nv)
plot(ecdf(tstat),verticals=TRUE,pch=NA,
     col='blue',xlim=c(-5,5),xlab='t',main='')
lines(ecdf(pred1),verticals=TRUE,pch=NA,col='black')
lines(ecdf(pred2),verticals=TRUE,pch=NA,col='red')
legend('topleft',legend=c('measured','n-1 d.o.f','n d.o.f'),
       lty=1,col=c('blue','black','red'))
dev.off()

cairo(file='KSsides.pdf',width=8,height=5)
set.seed(1)
pars(mfrow=c(3,4),mar=c(3,3,1,0))
ns <- 500
ni <- 1000
mA <- rep(10,3)
mB <- c(10,9,11)
ecdftitles <- c('A=B','A<B','A>B')
for (k in 1:3){
    Dm <- rep(Inf,ni)
    DM <- rep(-Inf,ni)
    D <- rep(Inf,ni)
    m <- rep(NA,ni)
    M <- rep(NA,ni)
    for (j in 1:ni){
        A <- rnorm(ns,mean=mA[k],sd=3)
        B <- rnorm(ns,mean=mB[k],sd=3)
        if (FALSE){
            cdfA <- ecdf(A)
            cdfB <- ecdf(B)
            for (i in 1:ns){
                a <- A[i]
                d <- cdfA(a) - cdfB(a)
                if (d<Dm[j]){
                    m[j] <- a
                    Dm[j] <- d
                }
                if (d>DM[j]){
                    M[j] <- a
                    DM[j] <- d
                }
                b <- B[i]
                d <- cdfA(b) - cdfB(b)
                if (d<Dm[j]){
                    m[j] <- b
                    Dm[j] <- d
                }
                if (d>DM[j]){
                    M[j] <- b
                    DM[j] <- d
                }
            }
        }
        if (TRUE){
            D[j] <- ks.test(A,B,alternative='two.sided')$statistic
            Dm[j] <- ks.test(A,B,alternative='less')$statistic
            DM[j] <- ks.test(A,B,alternative='greater')$statistic
        }
    }
    plot(ecdf(A),verticals=TRUE,pch=NA,col='blue',xlim=range(c(A,B)),main='')
    legend('topleft',legend=c('A','B'),lty=rep(1,2),col=c('red','blue'),bty='n')
    lines(ecdf(B),verticals=TRUE,pch=NA,col='red')
    hist(Dm,xlab=expression('D'^'-'),main='')
    rug(Dm)
    hist(DM,xlab=expression('D'^'+'),main='')
    rug(DM)
    hist(D,xlab='D',main='')
    rug(D)
}
dev.off()
