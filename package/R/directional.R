#' @title plot circular data
#' @description Plots directional data as ticks on a circle
#' @details Produces a circle with angles plotting in a clockwise
#'     direction from the top
#' @param angles scalar or vector
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @param tl tick length (value between 0 and 1)
#' @param ... optional arguments to be passed on to the generic
#'     \code{matlines} function
#' @examples
#' data(striations,package='geostats')
#' circle.plot(striations,degrees=TRUE)
#' @export
circle.plot <- function(angles,degrees=FALSE,tl=0.1,...){
    graphics::plot(x=c(-1,1),y=c(-1,1),type='n',axes=FALSE,
                   ann=FALSE,asp=1,bty='n')
    graphics::symbols(0,0,circles=1,add=TRUE,inches=FALSE)
    circle.markers(tl=tl)
    circle.ticks(angles,degrees=degrees,tl=tl,...)
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

circle.ticks <- function(angles,degrees=FALSE,tl=0.1,...){
    if (degrees) rads <- angles*pi/180
    else rads <- angles
    x1 <- sin(rads)
    y1 <- cos(rads)
    x2 <- x1*(1-tl)
    y2 <- y1*(1-tl)
    graphics::matlines(x=rbind(x1,x2),y=rbind(y1,y2),
                       lty=1,col='black',...)
}

#' @title add points to a circular plot
#' @description adds directional data as points on an existing circle
#'     plot
#' @details adds points to a circle with angles plotting in a
#'     clockwise direction from the top
#' @param angles scalar or vector
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @param ... optional arguments to be passed on to the generic
#'     \code{points} function
#' @examples
#' data(striations,package='geostats')
#' circle.plot(striations,degrees=TRUE)
#' md <- meanangle(striations,degrees=TRUE)
#' circle.points(md,pch=22,bg='black')
#' @export
circle.points <- function(angles,degrees=FALSE,...){
    if (degrees) rads <- angles*pi/180
    else rads <- angles
    graphics::points(x=sin(rads),y=cos(rads),...)
}

#' @title von Mises distribution
#' @description returns the probability density of a von Mises distribution
#' @details the von Mises distribution describes probability distributions
#' on a circle using the following density function:
#'
#' \eqn{\frac{\exp(\kappa\cos(x-\mu))}{2\pi I_0(\kappa)}}
#'
#' where \eqn{I_0(\kappa)} is a zero order Bessel function
#' 
#' @param angles scalar or vector
#' @param mu scalar containing the mean direction
#' @param kappa scalar containing the concentration parameter
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @return a scalar or vector of the same length as \code{angles}
#' @examples
#' plot(x=c(-1.2,1.2),y=c(-1.2,1.2),type='n',
#'      axes=FALSE,ann=FALSE,bty='n',asp=1)
#' a <- seq(from=-pi,to=pi,length.out=200)
#' d <- vonMises(angle=a,mu=pi/2,kappa=5)
#' symbols(x=0,y=0,circles=1,add=TRUE,inches=FALSE,xpd=NA,fg='grey50')
#' lines(x=x0+(1+d)*cos(a),y=(1+d)*sin(a),xpd=NA)
#' @export
vonMises <- function(angle,mu=0,kappa=1,degrees=FALSE){
    if (degrees) rad <- angles*pi/180
    num <- exp(kappa*cos(angle-mu))
    den <- 2*pi*besselI(kappa,nu=0)
    num/den
}

#' @title von Mises distribution
#' @description returns the probability density of a von Mises distribution
#' @details the von Mises distribution describes probability distributions
#' on a circle using the following density function:
#'
#' \eqn{\frac{\exp(\kappa\cos(x-\mu))}{2\pi I_0(\kappa)}}
#'
#' where \eqn{I_0(\kappa)} is a zero order Bessel function
#' 
#' @param angles scalar or vector
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @return the mean angle, either in radians (if
#'     \code{degrees=FALSE}), or in degrees.
#' @examples
#' data(striations,package='geostats')
#' circle.plot(angles=striations,degrees=TRUE)
#' circle.points(meanangle(striations,degrees=TRUE),pch=19)
#' @export
meanangle <- function(angles,degrees=FALSE){
    if (degrees) rad <- angles*pi/180
    out <- atan(sum(sin(rad))/sum(cos(rad)))
    if (degrees) out <- out*180/pi
    out
}

#' @title calculate \eqn{\bar{R}}
#' @description returns \eqn{\bar{R}}, a measure of directional concentration
#' @details Given \eqn{n} directional measurements \eqn{\theta_i},
#' 
#' \eqn{
#' \bar{R} =
#' \sqrt{\frac{\sum_{i=1}^{n}(\sin(\theta_i)^2 + cos(\theta_i)^2))}{n} }
#' }
#' 
#' @param angles scalar or vector
#' @param degrees \code{TRUE} for degrees, \code{FALSE} for radians
#' @return a value between 0 and 1
#' @examples
#' data(striations,package='geostats')
#' Rbar(angles=striations,degrees=TRUE)
#' @export
Rbar <- function(angles,degrees=FALSE){
    if (degrees) rad <- angles*pi/180
    sqrt(sum(sin(rad)^2 + cos(rad)^2))/length(angles)
}

#' @title \eqn{\bar{R}} to \eqn{\kappa} conversion
#' @description converts concentration parameter \eqn{\bar{R}} to
#'     \eqn{\kappa}
#' @details \eqn{\bar{R}} and \eqn{\kappa} are two types of
#'     concentration parameter that are commonly used in directional
#'     data analysis.  \eqn{\kappa} is one of the parameters of the
#'     parametric von Mises distribution, which is difficult to
#'     estimate from the data. \eqn{\bar{R}} is easier to calculate
#'     from data. \code{R2kappa} converts \eqn{\bar{R}} to
#'     \eqn{\bar{\kappa}} using a lookup table.
#' 
#' @param R a scalar or vector of values between 0 and 1
#' @return value(s) between 0 and \eqn{+\infty}
#' @examples
#' data(striations,package='geostats')
#' Rbar2kappa(Rbar(striations,degrees=TRUE))
#' @export
Rbar2kappa <- function(R){
    stats::spline(x=seq(from=0,to=1,by=0.01),
                  y=c(.00000,.02000,.04001,.06003,.08006,.10013,.12022,.14034,
                      .16051,.18073,.20101,.22134,.24175,.26223,.28279,.30344,
                      .32419,.34503,.36599,.38707,.40828,.42962,.45110,.47273,
                      .49453,.51649,.53863,.56097,.58350,.60625,.62922,.65242,
                      .67587,.69958,.72356,.74783,.77241,.79730,.82253,.84812,
                      .87408,.90043,.92720,.95440,.98207,1.01022,1.03889,1.06810,
                      1.09788,1.12828,1.15932,1.19105,1.22350,1.25672,1.29077,
                      1.32570,1.36156,1.39842,1.43635,1.47543,1.51574,1.55738,
                      1.60044,1.64506,1.69134,1.73945,1.78953,1.84177,1.89637,
                      1.95357,2.01363,2.07685,2.14359,2.21425,2.28930,2.36930,
                      2.45490,2.54686,2.64613,2.75382,2.87129,3.00020,3.14262,
                      3.30114,3.47901,3.68041,3.91072,4.17703,4.48876,4.85871,
                      5.3047,5.8522,6.5394,7.4257,8.6104,10.2716,12.7661,16.9266,
                      25.2522,50.2421,100),
                  xout=R)$y
}
