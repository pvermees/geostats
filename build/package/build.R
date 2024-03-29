setwd('~/Documents/Programming/R/geostats/build/package/')

catchments <- data.frame(
    lithology=c('basalt','granite','basalt','sandstone','shale',
                'basalt','shale','basalt','shale','shale','sandstone',
                'shale','sandstone','basalt','sandstone','granite',
                'basalt','shale','granite','basalt'),
    age=c('Cenozoic','Mesozoic','Mesozoic','Cenozoic','Precambrian',
          'Mesozoic','Palaeozoic','Palaeozoic','Mesozoic','Mesozoic',
          'Palaeozoic','Precambrian','Mesozoic','Cenozoic','Palaeozoic',
          'Palaeozoic','Mesozoic','Palaeozoic','Mesozoic','Cenozoic'),
    springs=c(0,1,0,6,3,1,2,0,1,2,4,2,4,1,4,0,0,1,0,1),
    pH=c(6.2,4.4,5.6,5.2,4.5,5.4,4.8,5.9,3.8,4.0,
         5.0,4.1,5.2,5.5,5.3,4.6,5.7,4.6,4.6,5.6),
    CaMg=c(0.35,11.00,6.00,1.80,2.30,0.59,8.40,2.90,5.90,2.10,
           1.20,2.10,1.10,1.60,0.90,1.70,3.40,0.53,2.20,7.70),
    vegetation=c(5.8,28.0,12,27,40,12,3.8,6.3,17,
                 16,95,94,92,88,88,70,92,72,74,84)
)
save(catchments,file="../../package/data/catchments.rda",version=2)

striations <- c(44,51,79,65,27,31,4,355,22,352,287,
                7,287,339,0,276,342,355,334,296,7,
                17,351,349,37,339,40,324,325,334)
save(striations,file="../../package/data/striations.rda",version=2)

pebbles <- c(32,33,34,37,40,41,42,43,45,53,
             210,212,214,217,220,221,222,223,226,230)
save(pebbles,file="../../package/data/pebbles.rda",version=2)

qdat <- read.csv('recentquakes.csv',header=TRUE)
dateAsString <- qdat$time
thedate <- strptime(dateAsString,format='%Y-%m-%d %H:%M:%S ')
earthquakes <- data.frame(
    year = as.numeric(format(thedate,'%Y')),
    month = as.numeric(format(thedate,'%m')),
    day = as.numeric(format(thedate,'%d')),
    hour = as.numeric(format(thedate,'%H')),
    minute = as.numeric(format(thedate,'%M')),
    second = as.numeric(format(thedate,'%S')),
    lat = qdat$latitude,
    lon = qdat$longitude,
    mag = qdat$mag
)
save(earthquakes,file='../../package/data/earthquakes.rda',compress='xz',version=2)

dqdat <- read.csv('declusteredquakes.csv',
                  header=TRUE)
declustered <- data.frame(
    year = dqdat$Year,
    month = dqdat$Month,
    day = dqdat$Day,
    hour = dqdat$Hour,
    minute = dqdat$Minute,
    second = dqdat$Second,
    lat = dqdat$Latitude,
    lon = dqdat$Longitude,
    mag = dqdat$Mwe
)
save(declustered,file='../../package/data/declustered.rda',
     version=2,compress='xz',compression_level=9)

dat <- read.csv('Finland.csv',header=TRUE)
Finland <- data.frame(
    area = dat$Lake_area,
    depth = dat$Depth_avg,
    lat = dat$Pour_lat,
    lon = dat$Pour_long,
    elevation = dat$Elevation
)
save(Finland,file='../../package/data/Finland.rda',
     version=2,compress='xz',compression_level=9)

set.seed(1)
DZ <- provenance::read.distributional('DZages.csv',check.names=FALSE)$x
nr <- max(unlist(lapply(DZ,'length')))
nc <- length(DZ)
DZages <- matrix('',nr,nc)
colnames(DZages) <- names(DZ)
done <- NULL
for (i in 1:nc){
    ng <- length(DZ[[i]])
    done <- c(done,DZ[[i]])
    nd <- length(done)
    duplicates <- duplicated(done)[(nd-ng+1):nd]
    DZij <- jitter(DZ[[i]])
    DZ[[i]][duplicates] <- round(DZij[duplicates],2)
    done[(nd-ng+1):nd] <- DZ[[i]]
    duplicates <- duplicated(done)[(nd-ng+1):nd]
    DZ[[i]][duplicates] <- round(DZij[duplicates],3)
    done[(nd-ng+1):nd] <- DZ[[i]]
    duplicates <- duplicated(done)[(nd-ng+1):nd]
    DZ[[i]][duplicates] <- round(DZij[duplicates],4)
    done[(nd-ng+1):nd] <- DZ[[i]]
    DZages[1:ng,i] <- DZ[[i]]
}
save(DZ,file="../../package/data/DZ.rda",version=2)
colnames(DZages) <- names(DZ)
write.csv(DZages,file='DZages.csv',row.names=FALSE)

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
        e <- MASS::mvrnorm(1,rep(0,2),E)
        RbSr[i] <- RbSr[i] + e[1]
        SrSr[i] <- SrSr[i] + e[2]
        errRbSr[i] <- sqrt(E[1,1])
        errSrSr[i] <- sqrt(E[2,2])
    }
    signif(rbind(RbSr,errRbSr,SrSr,errSrSr,rho),3)    
}

rbsr <- data.frame(t(RbSrGenerator(n=8)))
save(rbsr,file='../../package/data/rbsr.rda',version=2)

forams <- rbind(c(9,1,13,15,16,10,10),
                c(20,15,35,30,40,20,18))
colnames(forams) <- c('uvula','scitula','quinqueloba',
                      'pachyderma','incompta',
                      'glutinata','bulloides')
rownames(forams) <- c('A','B')
save(forams,file='../../package/data/forams.rda',version=2)

worldpop <- read.csv('population.csv')
save(worldpop,file='../../package/data/worldpop.rda',version=2)

raster2dat <- function(fname){
    mat <- dat <- tiff::readTIFF(fname)
    if (length(dim(dat))>2) mat <- dat <- dat[,,1]
    mat[dat<.5] <- 1
    mat[dat>=.5] <- 0
    t(apply(mat, 2, rev))
}

fractures <- raster2dat('fractures.tif')
save(fractures,file='../../package/data/fractures.rda',
     version=2,compress='bzip2',compression_level=9)

Corsica <- raster2dat('Corsica.tif')
save(Corsica,file='../../package/data/Corsica.rda',version=2)

Britain <- raster2dat('Britain.tif')
save(Britain,file='../../package/data/Britain.rda',version=2)

ACNK <- read.csv('ACNK.csv',header=TRUE,row.names=1,check.names=FALSE)
save(ACNK,file='../../package/data/ACNK.rda',version=2)

major <- read.csv('Major.csv',header=TRUE,row.names=1,check.names=FALSE)
save(major,file='../../package/data/major.rda',version=2)

FAM <- read.csv('FAM.csv',header=TRUE,check.names=FALSE)
save(FAM,file='../../package/data/FAM.rda',version=2)

dat <- read.delim('training.txt',header=TRUE,sep='\t',check.names=FALSE)
training <- na.omit(dat[,c(1,5:7,10:14)])
colnames(training) <- c('affinity','SiO2','TiO2','Al2O3','CaO',
                     'MgO','MnO','K2O','Na2O')
rownames(training) <- 1:nrow(training)
save(training,file='../../package/data/training.rda',compress='xz',version=2)

dat <- read.delim('test.txt',header=TRUE,sep='\t',check.names=FALSE)
test <- na.omit(dat[,c(1,5:7,10:14)])
colnames(test) <- c('affinity','SiO2','TiO2','Al2O3','CaO',
                    'MgO','MnO','K2O','Na2O')
rownames(test) <- 1:nrow(test)
save(test,file='../../package/data/test.rda',version=2)

palaeomag <- data.frame(decl=c(47.9,46.3,44.7,50.9,56.4,42.6,44.9,41.5,47.9,39.6),
                        incl=c(28.6,20.1,15.6,18.1,17.5,28.7,12.2,24.5,20.6,15.0))
save(palaeomag,file='../../package/data/palaeomag.rda',version=2)

fault <- data.frame(strike=c(311,319,316,319,324,312,314,319,320,314),
                    dip=c(38.3,36.5,34.2,34,35.5,32.9,40.7,37.5,40.5,43.2))
save(fault,file='../../package/data/fault.rda',version=2)

data('meuse',package='sp')
meuse <- meuse[c('x','y','elev','cadmium','lead','zinc')]
save(meuse,file='../../package/data/meuse.rda',version=2)

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
hills <- data.frame(X=X,Y=Y,Z=Z)
save(hills,file='../../package/data/hills.rda',version=2)
