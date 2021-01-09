#' @title
#' Linear regression of X,Y-variables with correlated errors
#'
#' @description
#' Implements the unified regression algorithm of York et al. (2004)
#' which, although based on least squares, yields results that are
#' consistent with maximum likelihood estimates of Titterington and
#' Halliday (1979).
#'
#' @details
#' Given n pairs of (approximately) collinear measurements \eqn{X_i}
#' and \eqn{Y_i} (for \eqn{1 \leq i \leq n}), their uncertainties
#' \eqn{s[X_i]} and \eqn{s[Y_i]}, and their covariances
#' cov[\eqn{X_i,Y_i}], the \code{york} function finds the best fitting
#' straight line using the least-squares algorithm of York et
#' al. (2004). This algorithm is modified from an earlier method
#' developed by York (1968) to be consistent with the maximum
#' likelihood approach of Titterington and Halliday (1979).
#'
#' @param dat a 4 or 5-column matrix with the X-values, the analytical
#'     uncertainties of the X-values, the Y-values, the analytical
#'     uncertainties of the Y-values, and (optionally) the correlation
#'     coefficients of the X- and Y-values.
#' @param alpha cutoff value for confidence intervals.
#' @param plot logical. If true, creates a scatter plot of the data
#'     with the best fit line shown on it.
#' @param fill the fill colour of the error ellipses. For additional
#'     plot options, use the \code{IsoplotR} package.
#' @param ... optional arguments for the scatter plot.
#' @return A two-element list of vectors containing:
#'
#'     \describe{
#'
#'     \item{coef}{the intercept and slope of the straight line fit}
#'
#'     \item{cov}{the covariance matrix of the coefficients}
#'
#'     }
#' 
#' @references
#' Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
#' parallel isochrons and the method of maximum likelihood. Chemical
#' Geology, 26(3), pp.183-195.
#'
#' York, Derek, et al., 2004. Unified equations for the slope,
#' intercept, and standard errors of the best straight line.  American
#' Journal of Physics 72.3, pp.367-375.
#'
#' @examples
#' data(rbsr,package='geostats')
#' fit <- york(rbsr)
#' @export
york <- function(dat,alpha=0.05,plot=TRUE,fill=NA,...){
    if (ncol(dat)==4) dat <- cbind(dat,0)
    colnames(dat) <- c('X','sX','Y','sY','rXY')
    ab <- stats::lm(dat[,'Y'] ~ dat[,'X'])$coefficients # initial guess
    a <- ab[1]
    b <- ab[2]
    if (any(is.na(ab)))
        stop('Cannot fit a straight line through these data')
    wX <- 1/dat[,'sX']^2
    wY <- 1/dat[,'sY']^2
    for (i in 1:50){ # 50 = maximum number of iterations
        bold <- b
        A <- sqrt(wX*wY)
        W <- wX*wY/(wX+b*b*wY-2*b*dat[,'rXY']*A)
        Xbar <- sum(W*dat[,'X'],na.rm=TRUE)/sum(W,na.rm=TRUE)
        Ybar <- sum(W*dat[,'Y'],na.rm=TRUE)/sum(W,na.rm=TRUE)
        U <- dat[,'X']-Xbar
        V <- dat[,'Y']-Ybar
        B <- W*(U/wY+b*V/wX-(b*U+V)*dat[,'rXY']/A)
        b <- sum(W*B*V,na.rm=TRUE)/sum(W*B*U,na.rm=TRUE)
        if ((bold/b-1)^2 < 1e-15) break # convergence reached
    }
    a <- Ybar-b*Xbar
    X <- Xbar + B
    xbar <- sum(W*X,na.rm=TRUE)/sum(W,na.rm=TRUE)
    u <- X-xbar
    sb <- sqrt(1/sum(W*u^2,na.rm=TRUE))
    sa <- sqrt(1/sum(W,na.rm=TRUE)+(xbar*sb)^2)
    out <- list()
    out$coef <- c(a,b)
    names(out$coef) <- c('intercept','slope')
    out$cov <- matrix(0,2,2)
    out$cov[1,1] <- sa^2
    out$cov[2,2] <- sb^2
    out$cov[1,2] <- -xbar*sb^2
    out$cov[2,1] <- -xbar*sb^2
    if (plot) scatterplot(dat,fit=out,fill=fill,...)
    invisible(out)
}

# get fitted X and X given a dataset x=cbind(X,sX,Y,sY,rXY),
# an intercept a and slope b. This function is useful
# for evaluating log-likelihoods of derived quantities
get.york.xy <- function(dat,a,b){
    wX <- 1/dat[,'sX']^2
    wY <- 1/dat[,'sY']^2
    A <- sqrt(wX*wY)
    W <- wX*wY/(wX+b*b*wY-2*b*dat[,'rXY']*A)
    Xbar <- sum(W*dat[,'X'],na.rm=TRUE)/sum(W,na.rm=TRUE)
    Ybar <- a + b*Xbar
    U <- dat[,'X']-Xbar
    V <- dat[,'Y']-Ybar
    B <- W*(U/wY+b*V/wX-(b*U+V)*dat[,'rXY']/A)
    out <- cbind(Xbar+B,Ybar+b*B)
    out
}

#' @title ellipse
#' @description Compute the x-y coordinates of an error ellipse.
#' @param mean two-element vector with the centre of the ellipse
#' @param cov the \code{2 x 2} covariance matrix of \code{x} and \code{y}
#' @param alpha confidence level of the ellipse
#' @param n the number of points at which the ellipse is evaluated
#' @return a two-column matrix of plot coordinates
#' @examples
#' X <- rnorm(100,mean=100,sd=1)
#' Y <- rnorm(100,mean=100,sd=1)
#' Z <- rnorm(100,mean=100,sd=5)
#' dat <- cbind(X/Z,Y/Z)
#' plot(dat)
#' ell <- ellipse(mean=colMeans(dat),cov=cov(dat))
#' polygon(ell)
#' @export
ellipse <- function(mean,cov,alpha=0.05,n=50){
    mu <- as.numeric(mean)
    cutoff <- stats::qchisq(1-alpha,2)
    e <- eigen(cov)
    a <- sqrt(cutoff*abs(e$values[1])) # major axis
    b <- sqrt(cutoff*abs(e$values[2])) # minor axis
    v <- e$vectors[,1] # largest eigenvector
    beta <- atan(v[2]/v[1]) # rotation angle of the ellipse
    theta <- seq(0, 2 * pi, length=n)
    out <- matrix(0,nrow=n,ncol=2)
    out[,1] <- mu[1] + a * cos(theta) * cos(beta) - b * sin(theta) * sin(beta)
    out[,2] <- mu[2] + a * cos(theta) * sin(beta) + b * sin(theta) * cos(beta)
    out
}

scatterplot <- function(xy,alpha=0.05,show.numbers=FALSE,
                        show.ellipses=1,levels=NA,clabel="",
                        fill=NA,stroke="black",fit='none',add=FALSE,
                        empty=FALSE, ci.col='gray80',line.col='black',
                        lwd=1,addcolourbar=TRUE,
                        bg,cex,xlim,ylim,xlab,ylab){
    ns <- nrow(xy)
    if (ncol(xy)==4) xy <- cbind(xy,rep(0,ns))
    colnames(xy) <- c('X','sX','Y','sY','rXY')
    if (missing(xlim)) xlim <- get.limits(xy[,'X'],xy[,'sX'])
    if (missing(ylim)) ylim <- get.limits(xy[,'Y'],xy[,'sY'])
    if (!add){
        if (missing(xlab)) xlab <- ''
        if (missing(ylab)) ylab <- ''
        graphics::plot(xlim,ylim,type='n',xlab=xlab,ylab=ylab,bty='n')
        if (empty) return()
    }
    if (!identical(fit,'none')){
        plot_isochron_line(fit,x=seq(xlim[1],xlim[2],length.out=100),
                           ci.col=ci.col,col=line.col,lwd=lwd)
    }
    graphics::box()
    nolevels <- all(is.na(levels))
    if (missing(cex)) cex <- 0.25
    for (i in 1:ns){
        if (!any(is.na(xy[i,]))){
            covmat <- cor2cov2(xy[i,'sX'],xy[i,'sY'],xy[i,'rXY'])
            ell <- ellipse(mean=xy[i,c('X','Y')],cov=covmat,alpha=alpha)
            graphics::polygon(ell,col=fill,border=stroke)
            if (show.numbers) graphics::text(xy[i,'X'],xy[i,'Y'],i)
            else graphics::points(xy[i,'X'],xy[i,'Y'],pch=19,cex=cex)
        }
    }
}

plot_isochron_line <- function(fit,x,ci.col='gray80',...){
    y <- fit$coef[1]+fit$coef[2]*x
    e <- sqrt(fit$cov[1,1] + 2*x*fit$cov[1,2] + fit$cov[2,2]*x^2)
    cix <- c(x,rev(x))
    ciy <- c(y+e,rev(y-e))
    graphics::polygon(cix,ciy,col=ci.col,border=NA)
    graphics::lines(x,y,...)
}

get.limits <- function(x,sx){
    minx <- min(x-3*sx,na.rm=TRUE)
    maxx <- max(x+3*sx,na.rm=TRUE)
    c(minx,maxx)
}

cor2cov2 <- function(sX,sY,rXY){
    covmat <- matrix(0,2,2)
    covmat[1,1] <- sX^2
    covmat[2,2] <- sY^2
    covmat[1,2] <- rXY*sX*sY
    covmat[2,1] <- covmat[1,2]
    covmat
}
