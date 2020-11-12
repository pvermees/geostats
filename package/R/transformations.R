#' @title logistic transformation
#' @description maps numbers from [0,1] to [\eqn{-\infty,+\infty}] and
#'     back
#' @param x a vector of real numbers (strictly positive if
#'     \code{inverse=FALSE})
#' @param inverse logical. If \code{inverse=FALSE}, returns
#'     \eqn{\ln\!\left[\frac{x}{1-x}\right]}; otherwise returns
#'     \eqn{\frac{\exp[x]}{\exp[x]+1}}.
#' @param ... optional arguments to the \code{log} function.
#' @return a vector with the same length of \code{x}
#' @examples
#' data(porosity,package='geostats')
#' lp <- logit(porosity,inverse=FALSE)
#' ld <- density(lp)
#' d <- logit(ld,inverse=TRUE)
#' plot(d)
#' @export
logit <- function(x,...){ UseMethod("logit",x) }
#' @rdname logit
#' @export
logit.default <- function(x,inverse=FALSE,...){
    if (inverse){
        out <- exp(x)/(exp(x)+1)
    } else if (all(x>0)){
        out <- log(x/(1-x),...)
    } else {
        stop('Cannot apply the logit transformation to negative values.')
    }
    out
}
#' @rdname logit
#' @export
logit.density <- function(x,inverse=TRUE,...){
    if (inverse | all(x$x>0)){
        out <- x
        out$x <- logit(x$x,inverse=inverse,...)
        dx <- diff(out$x)
        out$y <- x$y/c(dx,utils::tail(dx,1))
    } else {
        stop('Cannot apply the logit transformation to negative values.')
    }
    out
}

#' exponential transformation
#'
#' Map the input from [\eqn{-\infty,+\infty}] to [\eqn{0,\infty}] by
#' taking exponents
#' @param x an object of class \code{density}
#' @return an object of class \code{density}
#' @examples
#' data(clasts,package='geostats')
#' lc <- log(clasts)
#' ld <- density(lc)
#' d <- exp(ld)
#' plot(d)
#' @name exp
#' @export
exp.density <- function(x){
    out <- x
    out$x <- exp(x$x)
    dx <- diff(out$x)
    out$y <- x$y/c(dx,utils::tail(dx,1))
    out    
}

#' @title centred logratio transformation
#' @description maps compositional data from an n-dimensional simplex
#'     to an n-dimensional Euclidean space with Aitchison's centred
#'     logratio transformation
#' @param dat an \code{n x m} matrix
#' @param inverse if \code{TRUE}, applies the inverse clr tranformation
#' @return an \code{n x m} matrix
#' @examples
#' xyz <- rbind(c(0.03,99.88,0.09),
#'              c(70.54,25.95,3.51),
#'              c(72.14,26.54,1.32))
#' colnames(xyz) <- c('a','b','c')
#' rownames(xyz) <- 1:3
#' pc <- prcomp(clr(xyz))
#' biplot(pc)
#' @export
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

#' @title additive logratio transformation
#' @description maps compositional data from an n-dimensional simplex
#'     to an (n-1)-dimensional Euclidean space with Aitchison's
#'     additive logratio transformation
#' @param dat an \code{n x m} matrix
#' @param inverse if \code{TRUE}, applies the inverse alr
#'     tranformation
#' @return if \code{inverse=FALSE}, returns an \code{(n-1) x m} matrix
#'     of logratios; otherwise returns an \code{(n+1) x m} matrix of
#'     compositional data whose columns add up to 1.
#' @examples
#' xyz <- rbind(c(0.03,99.88,0.09),
#'              c(70.54,25.95,3.51),
#'              c(72.14,26.54,1.32))
#' colnames(xyz) <- c('a','b','c')
#' rownames(xyz) <- 1:3
#' pc <- prcomp(alr(xyz))
#' biplot(pc)
#' @export
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
