#' @title Principal Component Analysis of 2D data
#' @description produces a 4-panel summary plot for two dimensional
#'     PCA for didactical purposes
#' @param X a matrix with two columns
#' @examples
#' X <- rbind(c(-1,7),c(3,2),c(4,3))
#' colnames(X) <- c('a','b')
#' PCA2D(X)
#' @export
PCA2D <- function(X){
    pc <- stats::prcomp(X)
    # calculate the end-points of two lines marking the principal components (PC):
    CL <- matrix(NA,4,2) # initialise the matrix of coordinates
    CL[1:2,] <- matrix(1,2,1) %*% pc$center + diag(pc$sdev) %*% t(pc$rotation)
    CL[3:4,] <- matrix(1,2,1) %*% pc$center - diag(pc$sdev) %*% t(pc$rotation)
    # set up the 4-panel plot:
    p <- graphics::par(mfrow=c(2,2),mar=c(5,5,1,1),xpd=TRUE)
    # initialise the 1st panel:
    rx <- range(X[,'a'],CL[,1]) # range of x-values
    ry <- range(X[,'b'],CL[,2]) # range of y-values
    graphics::plot(rx,ry,type='n',asp=1,xlab='a',ylab='b')
    graphics::mtext('i',side=3,line=-1,adj=0.99)
    graphics::text(X,labels=1:3)
    # draw the line marking the 1st PC:
    graphics::lines(CL[c(1,3),])
    graphics::text(CL[3,1],CL[3,2],labels='PC1',pos=4)
    # draw the line marking the 2nd PC:
    graphics::lines(CL[c(2,4),])
    graphics::text(CL[2,1],CL[2,2],labels='PC2',pos=4)
    # add the centre point as a yellow square:
    graphics::points(t(pc$center),pch=22,bg='yellow')
    # initialise the 2nd panel:
    graphics::plot(range(pc$x),c(1,4),type='n',bty='n',
                   xaxt='n',yaxt='n',xlab='',ylab='')
    graphics::mtext('ii',side=3,line=-1,adj=0.99)
    graphics::Axis(side=1)
    # plot the 1st PC scores as a 1D configuration:
    graphics::lines(pc$x[,'PC1'],rep(2,3))
    graphics::points(pc$x[,'PC1'],rep(2,3))
    graphics::text(pc$x[,'PC1'],rep(2,3),labels=1:3,pos=c(1,1,3))
    graphics::text(min(pc$x[,'PC1']),2,labels='PC1',pos=2)
    # plot the 2nd PC scores as a 1D configuration:
    graphics::lines(pc$x[,'PC2'],rep(3,3))
    graphics::points(pc$x[,'PC2'],rep(3,3))
    graphics::text(pc$x[,'PC2'],rep(3,3),labels=1:3,pos=1)
    graphics::text(min(pc$x[,'PC2']),3,labels='PC2',pos=2)
    # plot both PCA scores and the loadings in the 3rd panel:
    stats::biplot(pc)
    graphics::mtext('iii',side=3,line=-1,adj=0.99)
    # plot the weights of the PCs in the 4th panel:
    w <- pc$sdev^2
    names(w) <- colnames(pc$x)
    graphics::barplot(w)
    graphics::mtext('iv',side=3,line=-1,adj=0.99)
    graphics::par(p)
}

#' @title Kolmogorov-Smirnov distance matrix
#' @description fills a square matrix with Kolmogorov-Smirnov
#'     statistics 
#' @param dat a list of numerical data vectors
#' @examples
#' data(DZ,package='geostats')
#' d <- ksdist(DZ)
#' plot(cmdscale(d))
#' @export
ksdist <- function(dat){
    snames <- names(dat)
    ns <- length(snames)
    out <- matrix(0,ns,ns)
    rownames(out) <- snames
    colnames(out) <- snames
    for (i in 1:ns){
        for (j in 1:ns){
            if (i!=j){
                out[i,j] <- stats::ks.test(dat[[i]],dat[[j]])$statistic
            }
        }
    }
    stats::as.dist(out)
}
