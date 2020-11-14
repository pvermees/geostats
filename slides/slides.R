setwd('/home/pvermees/Dropbox/teaching/Stats/slides')

pars <- function(mar=c(2.5,2.3,0.5,0.2),mgp=c(1.5,0.5,0),mfrow=c(1,1)){
    par(list(mar=mar,mgp=mgp,mfrow=mfrow))
}

cairo(file='binom10HT.pdf',width=4,height=4)
bpdf <- dbinom(0:3,size=3,p=0.5)
names(bpdf) <- 0:3
barplot(bpdf,col='NA')
dev.off()

cairo(file='poisquakes.pdf',width=4,height=4)
ppdf <- dpois(0:20,lambda=5.43)
names(ppdf) <- 0:20
barplot(ppdf,col='NA')
dev.off()

cairo(file='clasts.pdf',width=4,height=3)
pars()
clasts <- c(10,5,6,20)
names(clasts) <- c('granite','basalt','gneiss','quartzite')
barplot(clasts/sum(clasts),col='white')
dev.off()

cairo(file='pH.pdf',width=4,height=3)
pars()
pH <- c(6.2,4.4,5.6,5.2,4.5,5.4,4.8,5.9,3.9,3.8,
        5.1,4.1,5.1,5.5,5.1,4.6,5.7,4.6,4.6,5.6)
hist(pH,freq=FALSE,breaks=seq(from=3,to=7,by=0.5),xlim=c(3,7),main='')
rug(pH)
dev.off()
