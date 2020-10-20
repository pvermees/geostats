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
        out$y <- x$y/c(dx,tail(dx,1))
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
    out$y <- x$y/c(dx,tail(dx,1))
    out    
}
