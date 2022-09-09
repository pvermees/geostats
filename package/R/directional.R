#' @title plot circular data
#' @description Plots directional data as ticks on a circle, with
#'     angles plotting in a clockwise direction from the top.
#' @param a angle(s), scalar or vector
#' @param degrees logical. \code{TRUE} for degrees, \code{FALSE} for
#'     radians
#' @param tl tick length (value between 0 and 1)
#' @param ... optional arguments to be passed on to the generic
#'     \code{matlines} function
#' @return no return value
#' @examples
#' data(striations,package='geostats')
#' circle.plot(striations,degrees=TRUE)
#' @export
circle.plot <- function(a,degrees=FALSE,tl=0.1,...){
    graphics::plot(x=c(-1,1),y=c(-1,1),type='n',axes=FALSE,
                   ann=FALSE,asp=1,bty='n')
    graphics::symbols(0,0,circles=1,add=TRUE,inches=FALSE)
    circle.markers(tl=tl)
    circle.ticks(a,degrees=degrees,tl=tl,...)
}

circle.markers <- function(tl=0.1){
    graphics::lines(c(0,0),c(1,1+tl))
    graphics::lines(c(0,0),-c(1,1+tl))
    graphics::lines(c(1,1+tl),c(0,0))
    graphics::lines(-c(1,1+tl),c(0,0))
    graphics::text(0,1+tl,labels='0',pos=3,xpd=NA,offset=0.1)
    graphics::text(1+tl,0,labels='90',pos=4,xpd=NA,offset=0.1)
    graphics::text(0,-1-tl,labels='180',pos=1,xpd=NA,offset=0.15)
    graphics::text(-1-tl,0,labels='270',pos=2,xpd=NA,offset=0.1)
}

circle.ticks <- function(a,degrees=FALSE,tl=0.1,...){
    if (degrees) rads <- a*pi/180
    else rads <- a
    x1 <- sin(rads)
    y1 <- cos(rads)
    x2 <- x1*(1-tl)
    y2 <- y1*(1-tl)
    graphics::matlines(x=rbind(x1,x2),y=rbind(y1,y2),
                       lty=1,col='black',...)
}

#' @title add points to a circular plot
#' @description Adds directional data as points on an existing circle
#'     plot, with angles plotting in a clockwise direction from the
#'     top.
#' @param a angle(s), scalar or vector
#' @param degrees logical. \code{TRUE} for degrees, \code{FALSE} for
#'     radians
#' @param ... optional arguments to be passed on to the generic
#'     \code{points} function
#' @return no return value
#' @examples
#' data(striations,package='geostats')
#' circle.plot(striations,degrees=TRUE)
#' md <- meanangle(striations,degrees=TRUE)
#' circle.points(md,pch=22,bg='black',degrees=TRUE)
#' @export
circle.points <- function(a,degrees=FALSE,...){
    if (degrees) rads <- a*pi/180
    else rads <- a
    graphics::points(x=sin(rads),y=cos(rads),...)
}

#' @title von Mises distribution
#' @description Returns the probability density of a von Mises
#'     distribution, which describes probability distributions on a
#'     circle using the following density function:
#' 
#'     \eqn{\frac{\exp(\kappa\cos(x-\mu))}{2\pi I_0(\kappa)}}
#'
#' where \eqn{I_0(\kappa)} is a zero order Bessel function.
#' 
#' @param a angle(s), scalar or vector
#' @param mu scalar containing the mean direction
#' @param kappa scalar containing the concentration parameter
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @return a scalar or vector of the same length as \code{angles}
#' @examples
#' plot(x=c(-1,1.2),y=c(-1,1.2),type='n',
#'      axes=FALSE,ann=FALSE,bty='n',asp=1)
#' a <- seq(from=-pi,to=pi,length.out=200)
#' d <- vonMises(a=a,mu=pi/4,kappa=5)
#' symbols(x=0,y=0,circles=1,add=TRUE,inches=FALSE,xpd=NA,fg='grey50')
#' lines(x=(1+d)*cos(a),y=(1+d)*sin(a),xpd=NA)
#' @export
vonMises <- function(a,mu=0,kappa=1,degrees=FALSE){
    if (degrees) {
        a <- a*pi/180
        mu <- mu*pi/180
    }
    num <- exp(kappa*cos(a-mu))
    den <- 2*pi*besselI(kappa,nu=0)
    num/den
}

meanorientation <- function(trd,plg=0,option=1,degrees=FALSE){
    xyz <- angle2coord(trd,plg,option=option,degrees=degrees)
    n <- nrow(xyz)
    m <- ifelse(option==0,2,3)
    B <- matrix(0,m,m)
    for (i in 1:n){
        B <- B + t(xyz[i,1:m,drop=FALSE])%*%xyz[i,1:m,drop=FALSE]
    }
    E <- eigen(B)
    uvw <- E$vectors[,1]
    if (option==0){
        out <- coord2angle(uvw[1],uvw[2],option=option,degrees=degrees)
    } else {
        out <- coord2angle(x=sign(uvw[3])*uvw[1],
                           y=sign(uvw[3])*uvw[2],
                           z=abs(uvw[3]),option=option,degrees=degrees)
    }
    out
}

#' @title mean angle
#' @description Computes the vector mean of a collection of circular
#'     measurements.
#' @param trd trend angle, in degrees, between 0 and 360 (if
#'     \code{degrees=TRUE}) or between 0 and \eqn{2\pi} (if
#'     \code{degrees=FALSE}).
#' @param plg (optional) plunge angle, in degrees, between 0 and 90
#'     (if \code{degrees=TRUE}) or between 0 and \eqn{2\pi} (if
#'     \code{degrees=FALSE}).
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @param option scalar. If \code{option=0}, then \code{plg} is
#'     ignored and the measurements are considered to be circular; if
#'     \code{option=1}, then \code{trd} is the azimuth and \code{plg}
#'     is the dip; if \code{option=2}, then \code{trd} is the strike
#'     and \code{plg} is the dip; if \code{option=3} then \code{trd}
#'     is the longitude and \code{plg} is the latitude.
#' @param orientation logical. If \code{TRUE}, estimates the mean
#'     angle by eigen decomposition rather than by vector
#'     summation. This is the right thing to do for orientation data
#'     in which, for example, an angle of 45 degrees is equivalent to
#'     an angle of 225 degrees.
#' @return a scalar of 2-element vector with the mean orientation,
#'     either in radians (if \code{degrees=FALSE}), or in degrees.
#' @examples
#' data(striations,package='geostats')
#' meanangle(striations,degrees=TRUE)
#' @export
meanangle <- function(trd,plg=0,option=0,degrees=FALSE,orientation=FALSE){
    if (orientation){
        out <- meanorientation(trd=trd,plg=plg,option=option,degrees=degrees)
    } else {
        out <- vectorsum(trd=trd,plg=plg,option=option,degrees=degrees,Rbar=FALSE)
    }
}

#' @title calculate \eqn{\bar{R}}
#' @description Given \eqn{n} circular or spherical measurements, the
#'     length of their normalised vector sum (\eqn{\bar{R}}) serves as
#'     a measure of directional concentration.
#' @param trd trend angle, in degrees, between 0 and 360 (if
#'     \code{degrees=TRUE}) or between 0 and \eqn{2\pi} (if
#'     \code{degrees=FALSE}).
#' @param plg (optional) plunge angle, in degrees, between 0 and 90
#'     (if \code{degrees=TRUE}) or between 0 and \eqn{2\pi} (if
#'     \code{degrees=FALSE}).
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @param option scalar. If \code{option=0}, then \code{plg} is
#'     ignored and the measurements are considered to be circular; if
#'     \code{option=1}, then \code{trd} is the azimuth and \code{plg}
#'     is the dip; if \code{option=2}, then \code{trd} is the strike
#'     and \code{plg} is the dip; if \code{option=3} then \code{trd}
#'     is the longitude and \code{plg} is the latitude.
#' @return a value between 0 and 1
#' @examples
#' data(striations,package='geostats')
#' Rbar(striations,degrees=TRUE)
#' @export
Rbar <- function(trd,plg=0,option=0,degrees=FALSE){
    vectorsum(trd=trd,plg=plg,option=option,degrees=degrees,Rbar=TRUE)
}

# helper function for meanangle and Rbar
vectorsum <- function(trd,plg=0,option=0,degrees=FALSE,Rbar=TRUE){
    xyz <- angle2coord(trd,plg,option=option,degrees=degrees)
    x <- xyz[,'x']
    y <- xyz[,'y']
    z <- xyz[,'z']
    R <- sqrt(sum(x)^2+sum(y)^2+sum(z)^2)
    if (Rbar) return(R/length(trd))
    xbar <- sum(x)/R
    ybar <- sum(y)/R
    zbar <- sum(z)/R
    coord2angle(xbar,ybar,zbar,option=option,degrees=degrees)
}

angle2coord <- function(trd,plg=0,option=0,degrees=FALSE){
    if (degrees){
        trd <- trd*pi/180
        plg <- plg*pi/180
    }
    if (option==0){
        x <- cos(trd)
        y <- sin(trd)
        z <- 0
    } else if (option==1){
        az <- trd
        dip <- plg
        x <- cos(dip)*cos(az)
        y <- cos(dip)*sin(az)
        z <- sin(dip)
    } else if (option==2){
        strike <- trd
        dip <- plg
        x <- -cos(dip)*sin(strike)
        y <- cos(dip)*cos(strike)
        z <- sin(dip)
    } else if (option==3){
        lon <- trd
        lat <- plg
        x <- cos(lat)*sin(lon)
        y <- sin(lat)
        z <- -cos(lat)*cos(lon)
    } else {
        stop('Illegal option supplied to stereonet.point')
    }
    return(cbind(x=x,y=y,z=z))
}

coord2angle <- function(x,y,z=0,option=0,degrees=FALSE){
    if (option==0){
        out <- atan(y/x)
    } else if (option==1){
        out <- cbind(atan(y/x),
                     asin(z))
        out[x<0,1] <- out[x<0,1] - pi
    } else if (option==2){
        out <- cbind(atan(-x/y),
                     asin(z))
        out[y<0,1] <- out[y<0,1] - pi
    } else {
        out <- cbind(atan(-x/z),
                     asin(y))
        out[z>0,1] <- out[z>0,1] - pi
    }
    if (degrees) return(out*180/pi)
    else return(out)
}

#' @title \eqn{\bar{R}} to \eqn{\kappa} conversion
#' @description Converts the empirical concentration parameter
#'     \eqn{\bar{R}} to the von-Mises concentration parameter
#'     \eqn{\kappa}.
#' @details \eqn{\bar{R}} and \eqn{\kappa} are two types of
#'     concentration parameter that are commonly used in directional
#'     data analysis.  \eqn{\kappa} is one of the parameters of the
#'     parametric von Mises distribution, which is difficult to
#'     estimate from the data. \eqn{\bar{R}} is easier to calculate
#'     from data. \code{Rbar2kappa} converts \eqn{\bar{R}} to
#'     \eqn{\bar{\kappa}} using the following approximate empirical
#'     formula:
#'
#' \eqn{\kappa =
#'     \frac{\bar{R}(p+1-\bar{R}^2)}{1-\bar{R}^2}
#' }
#'
#' where \eqn{p} marks the number of parameters in the data space (1
#' for circle, 2 for a sphere).
#' @param R a scalar or vector of values between 0 and 1
#' @param p the number of parameters
#' @return value(s) between 0 and \eqn{+\infty}
#' @references Banerjee, A., et al. ``Clustering on the unit
#'     hypersphere using von Mises-Fisher distributions.''  Journal of
#'     Machine Learning Research 6.Sep (2005): 1345-1382.
#' @examples
#' data(striations,package='geostats')
#' Rbar2kappa(Rbar(striations,degrees=TRUE))
#' @export
Rbar2kappa <- function(R,p=1){
    R*(p+1-R^2)/(1-R^2)
}
