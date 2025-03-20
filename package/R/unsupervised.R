#' @title Principal Component Analysis of 2D data
#' @description Produces a 4-panel summary plot for two dimensional
#'     PCA for didactical purposes.
#' @param X a matrix with two columns
#' @param pos1 the positions of the text labels for the first
#'     principal component in the second panel
#' @param pos2 the positions of the text labels for the second
#'     principal component in the second panel
#' @examples
#' X <- cbind(a=c(0.0,1.7,1.9),b=c(0.9,0.0,0.6))
#' PCA2D(X)
#' @export
PCA2D <- function(X,
                  pos1=rep(c(1,3),ceiling(nrow(X)/2))[1:nrow(X)],
                  pos2=rep(c(1,1,3,3),ceiling(nrow(X)/4))[1:nrow(X)]){
    nr <- nrow(X)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    pc <- stats::princomp(X)
    # calculate the end-points of two lines marking the principal components (PC):
    CL <- matrix(NA,4,2) # initialise the matrix of coordinates
    CL[1:2,] <- matrix(1,2,1) %*% pc$center + diag(pc$sdev) %*% t(pc$loadings)
    CL[3:4,] <- matrix(1,2,1) %*% pc$center - diag(pc$sdev) %*% t(pc$loadings)
    # set up the 4-panel plot:
    graphics::par(mfrow=c(2,2),mar=c(5,5,1,1),xpd=TRUE)
    # initialise the 1st panel:
    rx <- range(X[,'a'],CL[,1]) # range of x-values
    ry <- range(X[,'b'],CL[,2]) # range of y-values
    graphics::plot(rx,ry,type='n',asp=1,xlab='a',ylab='b')
    graphics::legend('topright',legend='i',bty='n')
    graphics::text(X,labels=1:nr)
    # draw the line marking the 1st PC:
    graphics::lines(CL[c(1,3),])
    graphics::text(CL[3,1],CL[3,2],labels='PC1',pos=4)
    # draw the line marking the 2nd PC:
    graphics::lines(CL[c(2,4),])
    graphics::text(CL[2,1],CL[2,2],labels='PC2',pos=4)
    # add the centre point as a yellow square:
    graphics::points(t(pc$center),pch=22,bg='yellow')
    # initialise the 2nd panel:
    graphics::plot(range(pc$scores),c(1,4),type='n',bty='n',
                   xaxt='n',yaxt='n',xlab='',ylab='')
    graphics::legend('topright',legend='ii',bty='n')
    graphics::Axis(side=1)
    # plot the 1st PC scores as a 1D configuration:
    graphics::lines(pc$scores[,1],rep(2,nr))
    graphics::points(pc$scores[,1],rep(2,nr))
    graphics::text(pc$scores[,1],rep(2,nr),labels=1:nr,pos=pos1)
    graphics::text(min(pc$scores[,1]),2,labels='PC1',pos=2)
    # plot the 2nd PC scores as a 1D configuration:
    graphics::lines(pc$scores[,2],rep(3,nr))
    graphics::points(pc$scores[,2],rep(3,nr))
    graphics::text(pc$scores[,2],rep(3,nr),labels=1:nr,pos=pos2)
    graphics::text(min(pc$scores[,2]),3,labels='PC2',pos=2)
    # plot both PCA scores and the loadings in the 3rd panel:
    stats::biplot(pc)
    graphics::legend('topright',legend='iii',bty='n')
    # plot the weights of the PCs in the 4th panel:
    w <- pc$sdev^2
    names(w) <- colnames(pc$scores)
    graphics::barplot(w)
    graphics::legend('topright',legend='iv',bty='n')
}

#' @title Kolmogorov-Smirnov distance matrix
#' @description Given a list of numerical vectors, fills a square
#'     matrix with Kolmogorov-Smirnov statistics.
#' @param dat a list of numerical data vectors
#' @return an object of class \code{dist}
#' @examples
#' data(DZ,package='geostats')
#' d <- ksdist(DZ)
#' mds <- cmdscale(d)
#' plot(mds,type='n')
#' text(mds,labels=names(DZ))
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
