#' @title get the mode of a dataset
#' @description Computes the most frequently occuring value in a
#'     sampling distribution.
#' @param x a vector
#' @param categorical logical. If \code{TRUE}, returns the most
#'     frequently occuring value for categorical variables. If
#'     \code{FALSE}, returns the value corresponding to the maximimum
#'     kernel density for continuous variables
#' @return a scalar
#' @examples
#' data(islands,package='geostats')
#' m1 <- Mode(islands$area,categorical=TRUE)
#' 
#' m2 <- 1:50
#' for (i in m2){
#'    m2[i] <- Mode(rnorm(100),categorical=FALSE)
#' }
#' hist(m2)
#' @export
Mode <- function(x,categorical=FALSE){
    if (categorical){
        uni <- unique(x)
        out <- uni[which.max(tabulate(match(x, uni)))]
    } else {
        dens <- stats::density(x)
        out <- dens$x[which.max(dens$y)]
    }
    out
}

#' @title calculate the skewness of a dataset
#' @description Compute the third moment of a sampling distribution.
#' @param x a vector
#' @return a scalar
#' @examples
#' data(islands,package='geostats')
#' skew(islands$vegetation)
#' @export
skew <- function(x){
    mean((x-mean(x))^3)/(length(x)*stats::sd(x)^3)
}

#' @title generate bivariate random data
#' @description Returns bivariate datasets from four synthetic
#'     distributions that have the shape of a circle, arrow, square
#'     and ellipse.
#' @param pop an integer from 1 to 4 marking the population of choice:
#' 1 = circle, 2 = arrow, 3 = solid square, 4 = ellipse.
#' @param n the number of random draws to be drawn from population \code{pop}
#' @return a \code{[2xn]} matrix of random numbers
#' @examples
#' p <- par(mfrow=c(1,4))
#' for (i in 1:4){
#'    plot(randy(pop=i))
#' }
#' par(p)
#' @export
randy <- function(pop=1,n=250){
    if (pop == 1){
        mu <- c(0.5,0.5)
        angle <- seq(0,2*pi,length.out=50)
        radius <- 0.5
        rangle <- stats::runif(n,min=0,max=2*pi)
        rx <- mu[1] + radius*cos(rangle)
        ry <- mu[2] + radius*sin(rangle)
    } else if (pop == 2){ # the arrow is made of three segments
        m <- 0
        M <- 1
        psegment1 <- sqrt(2)*(M-m)/(sqrt(2)*(M-m)+0.5) # the tail
        psegment2 <- 0.25/(sqrt(2)*(M-m)+0.5) # one half of the head
        cutoff1 <- psegment1
        cutoff2 <- psegment1 + psegment2
        rnum <- stats::runif(n)
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
        rx <- stats::runif(n,min=min(x),max=max(x))
        ry <- stats::runif(n,min=min(y),max=max(y))
    } else if (pop == 4){
        mu <- c(0.5,0.5)
        angle <- seq(0,2*pi,length.out=50)
        theta <- -pi/4
        a <- 0.2
        b <- 0.66
        rangle <- stats::runif(n,min=0,max=2*pi)
        rx <- mu[1] + a*cos(rangle)*cos(theta) - b*sin(rangle)*sin(theta)
        ry <- mu[2] + a*cos(rangle)*sin(theta) + b*sin(rangle)*cos(theta)
    }
    cbind(rx,ry)
}

hasClass <- function(obj,...){
    any(class(obj)%in%unlist(list(...)))
}
