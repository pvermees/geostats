# returns a table with the number of earthquakes in the USGS
# catalog ('qdat') whose magnitude exceeds 'mag'
countQuakes <- function(qdat,mag){
    largequakes <- (qdat$Magnitude > mag)
    dateAsString <- qdat$DateTime[largequakes]
    thedate <- strptime(dateAsString,format='%Y/%m/%d')
    yearAsString <- format(thedate,'%Y')
    yearAsNum <- as.numeric(yearAsString)
    nquakes <- table(yearAsNum)
    return(nquakes)
}

countFloods <- function(dat,flow=350,from=1900,to=2000){
    isflood <- dat[,'flow']>flow
    selection <- dat[,'year']>=from & dat[,'year']<=to
    years <- unique(dat[selection,'year'])
    ndecades <- length(years)
    out <- rep(0,ny)
    for (i in 1:ny){
        theyear <- years[i]
        thisyear <- which(dat[,'year'] %in% theyear)
        blocks <- rle(isflood[thisyear])
        out[i] <- sum(blocks$values)
    }
    names(out) <- years
    out
}

sierpinski <- function(level){
    IterateCarpet <- function(A){
        B <- cbind(A,A,A)
        C <- cbind(A,0*A,A)
        D <- rbind(B,C,B)
        return(D)
    }
    S <- matrix(1,1,1)
    for (i in 1:level)
        S <- IterateCarpet(S)
    return(S)
}

# pop = an integer from 1 to 4 marking the population of choice
# n = the number of random draws to be drawn from population pop
# returns a [2xn] matrix of random numbers
randy <- function(pop=1,n=250){
    if (pop == 1){
        mu <- c(0.5,0.5)
        angle <- seq(0,2*pi,length.out=50)
        radius <- 0.5
        rangle <- runif(n,min=0,max=2*pi)
        rx <- mu[1] + radius*cos(rangle)
        ry <- mu[2] + radius*sin(rangle)
    } else if (pop == 2){ # the arrow is made of three segments
        m <- 0
        M <- 1
        psegment1 <- sqrt(2)*(M-m)/(sqrt(2)*(M-m)+0.5) # the tail
        psegment2 <- 0.25/(sqrt(2)*(M-m)+0.5) # one half of the head
        cutoff1 <- psegment1
        cutoff2 <- psegment1 + psegment2
        rnum <- runif(n)
        fromtail <- (rnum < cutoff1)
        fromhead1 <- (rnum >= cutoff1 & rnum < cutoff2)
        fromhead2 <- (rnum > cutoff2)
        rx <- rep(0,n)
        ry <- rep(0,n)
        rx[fromtail] <-   m+(M-m)*rnum[fromtail]/cutoff1
        ry[fromtail] <-   m+(M-m)*rnum[fromtail]/cutoff1
        rx[fromhead1] <-  M
        ry[fromhead1] <-  0.5+0.5*(rnum[fromhead1]-cutoff1)/(cutoff2-cutoff1)
        rx[fromhead2] <-  0.5+0.5*(rnum[fromhead2]-cutoff2)/(1-cutoff2)
        ry[fromhead2] <-  M
    } else if (pop == 3){
        x <- c(0,1,1,0,0)
        y <- c(0,0,1,1,0)
        rx <- runif(n,min=min(x),max=max(x))
        ry <- runif(n,min=min(y),max=max(y))
    } else if (pop == 4){
        mu <- c(0.5,0.5)
        angle <- seq(0,2*pi,length.out=50)
        theta <- -pi/4
        a <- 0.2
        b <- 0.66
        rangle <- runif(n,min=0,max=2*pi)
        rx <- mu[1] + a*cos(rangle)*cos(theta) - b*sin(rangle)*sin(theta)
        ry <- mu[2] + a*cos(rangle)*sin(theta) + b*sin(rangle)*cos(theta)
    }
    cbind(rx,ry)
}

# add a 'confidence polygon' to an existing ternary diagram
ternary.polygon <- function(LL,UL,...){
    min.a <- LL[1]
    max.a <- UL[1]
    min.b <- LL[2]
    max.b <- UL[2]
    min.c <- LL[3]
    max.c <- UL[3]
    corners <- matrix(c(max.a,min.b,1-max.a-min.b,
                        max.a,1-max.a-min.c,min.c,
                        1-max.b-min.c,max.b,min.c,
                        min.a,max.b,1-min.a-max.b,
                        min.a,1-min.a-max.c,max.c,
                        1-min.b-max.c,min.b,max.c,
                        max.a,min.b,1-max.a-min.b),
                      byrow=TRUE,ncol=3)
    polygon <- provenance:::xyz2xy(corners)
    i <- chull(polygon)
    lines(polygon[c(i,i[1]),],...)
}

# 4-panel summary plot for 2D principal component analysis
# X = a 2-column table of numerical values
PCA2D <- function(X){
    pc <- prcomp(X)
    # calculate the end-points of two lines marking the principal components (PC):
    CL <- matrix(NA,4,2) # initialise the matrix of coordinates
    CL[1:2,] <- matrix(1,2,1) %*% pc$center + diag(pc$sdev) %*% t(pc$rotation)
    CL[3:4,] <- matrix(1,2,1) %*% pc$center - diag(pc$sdev) %*% t(pc$rotation)
    # set up the 4-panel plot:
    par(mfrow=c(2,2),mar=c(5,5,1,1),xpd=TRUE)
    # initialise the 1st panel:
    rx <- range(X[,'a'],CL[,1]) # range of x-values
    ry <- range(X[,'b'],CL[,2]) # range of y-values
    plot(rx,ry,type='n',asp=1,xlab='a',ylab='b')
    mtext('i',side=3,line=-1,adj=0.99)
    text(X,labels=1:3)
    # draw the line marking the 1st PC:
    lines(CL[c(1,3),])
    text(CL[3,1],CL[3,2],labels='PC1',pos=4)
    # draw the line marking the 2nd PC:
    lines(CL[c(2,4),])
    text(CL[2,1],CL[2,2],labels='PC2',pos=4)
    # add the centre point as a yellow square:
    points(t(pc$center),pch=22,bg='yellow')
    # initialise the 2nd panel:
    plot(range(pc$x),c(1,4),type='n',bty='n',
         xaxt='n',yaxt='n',xlab='',ylab='')
    mtext('ii',side=3,line=-1,adj=0.99)
    Axis(side=1)
    # plot the 1st PC scores as a 1D configuration:
    lines(pc$x[,'PC1'],rep(2,3))
    points(pc$x[,'PC1'],rep(2,3))
    text(pc$x[,'PC1'],rep(2,3),labels=1:3,pos=c(1,1,3))
    text(min(pc$x[,'PC1']),2,labels='PC1',pos=2)
    # plot the 2nd PC scores as a 1D configuration:
    lines(pc$x[,'PC2'],rep(3,3))
    points(pc$x[,'PC2'],rep(3,3))
    text(pc$x[,'PC2'],rep(3,3),labels=1:3,pos=1)
    text(min(pc$x[,'PC2']),3,labels='PC2',pos=2)
    # plot both PCA scores and the loadings in the 3rd panel:
    biplot(pc)
    mtext('iii',side=3,line=-1,adj=0.99)
    # plot the weights of the PCs in the 4th panel:
    w <- pc$sdev^2
    names(w) <- colnames(pc$x)
    barplot(w)
    mtext('iv',side=3,line=-1,adj=0.99)
}

pebbles <- c(44,51,79,65,27,31,4,355,22,352,287,
             7,287,339,0,276,342,355,334,296,7,
             17,351,349,37,339,40,324,325,334)

plot.circ <- function(angles,degrees=FALSE,...){
    buffer <- 0.1
    plot(x=c(-1-buffer,1+buffer),y=c(-1-buffer,1+buffer),
         type='n',axes=FALSE,ann=FALSE,asp=1)
    symbols(0,0,circles=1,add=TRUE,inches=FALSE)
    points.circ(angles,degrees=degrees,...)
    tl <- 0.05 # tick length
    lines(c(0,0),c(1,1+tl))
    lines(c(0,0),-c(1,1+tl))
    lines(c(1,1+tl),c(0,0))
    lines(-c(1,1+tl),c(0,0))
    text(0,1+tl,labels='0',pos=3)
    text(1+tl,0,labels='90',pos=4)
    text(0,-1-tl,labels='180',pos=1)
    text(-1-tl,0,labels='270',pos=2)
}

points.circ <- function(angles,degrees=FALSE,...){
    if (degrees) rads <- angles*pi/180
    else rads <- angles
    points(sin(rads),cos(rads),...)    
}

