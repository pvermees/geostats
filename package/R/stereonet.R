#' @title stereonet
#' @description Plots directional data on a Schmidt equal area or
#'     Wulff equal angle stereonet.
#' @details The Schmidt equal area polar Lambert projection is used to
#'     plot structural data. Given a set of azimuth (\eqn{a}) and dip
#'     (\eqn{d}) measurements, it plots the data at positions
#'
#' \eqn{x = \sqrt{2}\sin(\pi\frac{90-d}{360})
#'                  \cos(\pi\frac{90-a}{180})}
#'
#' and
#'
#' \eqn{y = \sqrt{2}\sin(\pi\frac{90-d}{360})
#'                  \sin(\pi\frac{90-a}{180})}
#'
#' The Wulff equal angle polar Lambert projection is used to plot
#'     crystallographic data. Given a set of azimuth (\eqn{a}) and dip
#'     (\eqn{d}) measurements, it plots the data at positions
#'
#' \eqn{x = \tan(\pi\frac{90-d}{360})
#'          \cos(\pi\frac{90-a}{180})}
#'
#' and
#'
#' \eqn{y = \tan(\pi\frac{90-d}{360})
#'          \sin(\pi\frac{90-a}{180})}
#' 
#' @param trd trend angle, in degrees, between 0 and 360 (if
#'     \code{degrees=TRUE}) or between 0 and \eqn{2\pi} (if
#'     \code{degrees=FALSE}).
#' @param plg plunge angle, in degrees, between 0 and 90 (if
#'     \code{degrees=TRUE}) or between 0 and \eqn{2\pi} (if
#'     \code{degrees=FALSE}).
#' @param option scalar. If \code{option=1), then \code{trd} is the
#'     azimuth and \code{plg} is the dip; if \code{option=2), then
#'     \code{trd} is the strike and \code{plg} is the dip; if
#'     \code{option=3), then \code{trd} is the longitude and
#'     \code{plg} is the latitude.
#' @param add logical. If \code{TRUE}, adds to an existing stereonet.
#' @param wulff logical. If \code{FALSE}, produces a Schmidt net.
#' @param degrees logical. If \code{FALSE}, assumes that
#'     \code{azimuth} and \code{dip} are in radians.
#' @param show.grid logical. If \code{TRUE}, decorates the plot with
#'     a grid of great and small circles.
#' @param tl tick length for the N, E, S, W markers (value between 0
#'     and 1).  Set to 0 to omit the markers.
#' @param ... optional arguments to be passed on to the generic
#'     \code{points} function
#' @author based on the Matlab script by Gerry Middleton
#' @examples
#' stereonet(azimuth=c(120,80),dip=c(10,30),degrees=TRUE)
#' @export
stereonet <- function(trd,plg,coneAngle=0,option=1,wulff=TRUE,add=FALSE,
                      degrees=FALSE,show.grid=TRUE,tl=0.05,...){
    if (!add){
        graphics::plot(x=c(-1-tl,1+tl),y=c(-1-tl,1+tl),type='n',asp=1,
                       bty='n',axes=FALSE,ann=FALSE)
        graphics::symbols(x=0,y=0,circles=1,add=TRUE,inches=FALSE)
        if (tl>0) circle.markers(tl=tl)
        if (show.grid){
            graphics::lines(x=c(-1,1),y=c(0,0))
            graphics::lines(x=c(0,0),y=c(-1,1))
            N <- 50
            psi = seq(from=0,to=pi,by=pi/N)
            for (i in 1:8){                          # plot great circles
                rdip <- i*(pi/18)                    # at 10 deg intervals
                radip <- atan(tan(rdip)*sin(psi))
                if (wulff) rproj <- tan((pi/2 - radip)/2)
                else rproj <- sqrt(2)*sin((pi/2 - radip)/2)
                x1 <- rproj * sin(psi)
                x2 <- rproj * (-sin(psi))
                y <- rproj * cos(psi)
                graphics::lines(x1,y,lty=3)
                graphics::lines(x2,y,lty=3)
            }
            for (i in 1:8){                          # plot small circles
                alpha <- i*(pi/18)
                xlim <- sin(alpha)
                ylim <- cos(alpha)
                x <- seq(from=-xlim,to=xlim,length.out=50)
                d <- 1/cos(alpha)
                rd <- d*sin(alpha)
                y0 <- sqrt(rd^2 - x^2)
                y1 <- d - y0
                y2 <- - d + y0
                graphics::lines(x,y1,lty=3)
                graphics::lines(x,y2,lty=3)
            }
        }
    }
    if (!missing(trd) & !missing(plg)){
        if (degrees){
            trd <- trd*pi/180
            plg <- plg*pi/180
        }
        if (option==1){
            stereonet.line(trd,plg,wulff=wulff,...)
        } else if (option==2){
            stereonet.plane(trd,plg,wulff=wulff,...)
        } else if (option==3){
            stereonet.point(trd,plg,wulff=wulff,...)
        } else if (option==4){
            paths <- SmallCircle(trda=trd,plga=plg,coneAngle=coneAngle,wulff=wulff)
            lines(paths$path1,...)
            lines(paths$path2,...)
        } else {
            stop('Illegal option in stereonet function.')
        }
    }
}

stereonet.line <- function(trd,plg,wulff=TRUE,...){
    theta <- (pi/2)-trd
    if (wulff) rho <- tan((pi/4)-(plg/2))
    else rho <- sqrt(2)*sin((pi/2 - plg)/2)
    xp <- rho*cos(theta)
    yp <- rho*sin(theta)
    graphics::points(xp,yp,...)
}
stereonet.plane <- function(trd,plg,wulff=TRUE,...){
    for (i in 1:length(trd)){
        ad <- pole(trd=trd[i],plg=plg[i],option=2)
        stereonet.line(trd=ad[1],plg=ad[2],wulff=wulff,...)
        xy <- GreatCircle(strike=trd[i],dip=plg[i],wulff=wulff)
        lines(xy,...)
    }
}

#' returns the pole to a plane or the plane which correspond to a pole
#' based on a MATLAB script written by Nestor Cardozo for the book
#' Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher,
#' 2011.  trd and plg are in radians
pole <- function(trd,plg,option=1){
    out <- c(trd,plg)
    east <- pi/2
    if (option == 1){ # azimuth to strike
        if (plg >= 0.0){
            out[2] <- east - plg
            dipaz <- trd - pi
        } else {
            out[2] <- east + plg
            dipaz <- trd
        }
        out[1] <- ZeroTwoPi(dipaz - east)
    } else if (option == 2){ # strike to azimuth
        cned <- SphToCart(trd,plg,option=option)
        out <- CartToSph(cn=cned[1],ce=cned[2],cd=cned[3]);
    } else {
        warning('Invalid value for option in pole function')
    }
    out
}

# ZeroTwoPi constrains azimuth to lie between 0 and 2*pi radians 
# based on a MATLAB script written by Nestor Cardozo for the book
# Structural Geology Algorithms by Allmendinger, Cardozo & Fisher, 2011.
ZeroTwoPi <- function(a){
    twopi <- 2*pi
    if (a < 0){
        out <- a + twopi
    } else if (a >= twopi){
        out <- a - twopi;
    } else {
        out <- a
    }
    out
}

# SphToCart converts from spherical to cartesian coordinates 
# based on a MATLAB script written by Nestor Cardozo for the book
# Structural Geology Algorithms by Allmendinger, Cardozo & Fisher, 2011.
SphToCart <- function(trd,plg,option=1){
    if (option == 1){ # trend and plunge of a line
        cd <- sin(plg);
        ce <- cos(plg) * sin(trd);
        cn <- cos(plg) * cos(trd); 
    } else if (option == 2){ # strike and dip of a plane in the right hand rule
        cd <- cos(plg);
        ce <- -sin(plg) * cos(trd);
        cn <- sin(plg) * sin(trd);
    } else {
        stop('Invalid value for option in SphToCart')
    }
    c(cn,ce,cd)
}

#' Convert Cartesian to spherical coordinates
#' based on a MATLAB script written by Nestor Cardozo for the book
#' Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011.
#' returns the trend (trd) and plunge (plg) of a line for
#' input north (cn), east (ce), and down (cd) direction cosines
CartToSph <- function(cn,ce,cd){
    plg <- asin(cd)
    # If north direction cosine is zero, trend is east or west
    # Choose which one by the sign of the east direction cosine
    if (cn == 0){
        if (ce < 0)
            trd <- 3/2*pi # trend is west
        else
            trd <- pi/2   # trend is east
    } else {
        trd <- atan(ce/cn)
        if (cn < 0){
            trd <- trd + pi
        }
        # Make sure trd is between 0 and 2*pi
        trd <- ZeroTwoPi(trd)
    }
    c(trd,plg)
}

#' computes the great circle path of a plane in an equal angle or
#' equal area stereonet of unit radius.  Basedon the MATLAB script
#' written by Nestor Cardozo for the book Structural Geology
#' Algorithms by Allmendinger, Cardozo, & Fisher, 2011.
GreatCircle <- function(strike,dip,wulff=TRUE){
    tpa <- pole(strike,dip,option=2)
    trd <- strike
    plg <- 0
    rot <- (0:180)*pi/180
    out <- matrix(0,nrow=180,ncol=2)
    for (i in 1:180){
        if (rot[i] == pi){
            rot[i] <- rot[i]*0.9999
        }
        rtp <- Rotate(tpa[1],tpa[2],rot[i],trd,plg,'a')
        out[i,] <- StCoordLine(rtp[1],rtp[2],wulff=wulff)
    }
    out
}

#' Computes the coordinates of a line on a stereonet.  Based on a
#' MATLAB script written by Nestor Cardozo for the book Structural
#' Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011.
StCoordLine <- function(trd,plg,wulff=TRUE){
    if (plg < 0){
        trd <- ZeroTwoPi(trd+pi)
        plg <- -plg
    }
    piS4 <- pi/4
    s2 <- sqrt(2)
    plgS2 <- plg/2
    if (wulff){
        xp <- tan(piS4 - plgS2)*sin(trd)
        yp <- tan(piS4 - plgS2)*cos(trd)
    } else {
        xp <- s2*sin(piS4 - plgS2)*sin(trd)
        yp <- s2*sin(piS4 - plgS2)*cos(trd)
    }
    cbind(xp,yp)
}

#' Rotate a line by performing a coordinate transformation on
#' vectors. The algorithm was originally written by Randall A. Marrett
#' and implemented in MATLAB by Nestor Cardozo for the book Structural
#' Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011.
#' raz = trend of rotation axis
#' rdip = plunge of rotation axis
#' rot = magnitude of rotation
#' trd = trend of the vector to be rotated
#' plg = plunge of the vector to be rotated
#' ans0 = a character indicating whether the line to be rotated
#' is an axis (ans0 = 'a') or a vector (ans0 = 'v')
Rotate <- function(raz,rdip,rot,trd,plg,ans0){
    a <- matrix(0,3,3)
    pole <- rep(0,3)
    plotr <- rep(0,3)
    temp <- rep(0,3)
    p <- SphToCart(raz,rdip,option=1)
    x <- 1 - cos(rot)
    sinRot <- sin(rot)
    cosRot <- cos(rot)
    a[1,1] <- cosRot + p[1]*p[1]*x
    a[1,2] <- -p[3]*sinRot + p[1]*p[2]*x
    a[1,3] <- p[2]*sinRot + p[1]*p[3]*x
    a[2,1] <- p[3]*sinRot + p[2]*p[1]*x
    a[2,2] <- cosRot + p[2]*p[2]*x
    a[2,3] <- -p[1]*sinRot + p[2]*p[3]*x
    a[3,1] <- -p[2]*sinRot + p[3]*p[1]*x
    a[3,2] <- p[1]*sinRot + p[3]*p[2]*x
    a[3,3] <- cosRot + p[3]*p[3]*x
    temp <- SphToCart(trd,plg,option=1)
    for (i in 1:3){
        plotr[i] = 0
        for (j in 1:3){
            plotr[i] = a[i,j]*temp[j] + plotr[i]
        }
    }
    if (plotr[3] < 0 && ans0 == 'a'){
        plotr[1] <- -plotr[1]
        plotr[2] <- -plotr[2]
        plotr[3] <- -plotr[3]
    }
    CartToSph(plotr[1],plotr[2],plotr[3])
}

#' Compute the paths of a small circle defined by its axis and cone
#' angle, for an equal angle or equal area stereonet of unit radius.
#' Based on a MATLAB script written by Nestor Cardozo for the book
#' Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher,
#' 2011.
SmallCircle <- function(trda,plga,coneAngle,wulff=TRUE){
    if (plga < coneAngle){
        if (plga == pi/2){
            plga = plga * 0.9999
        }
        angle <- acos(cos(coneAngle)/cos(plga))
        trd <- ZeroTwoPi(trda+angle)
        plg <- 0
    } else {
        trd <- trda
        plg <- plga - coneAngle
    }
    rot <- (0:360)*pi/180
    path1 <- NULL
    path2 <- NULL
    np1 <- 0
    np2 <- 0
    for (i in 1:360){
        rtp = Rotate(trda,plga,rot[i],trd,plg,'v')
        if (rtp[2] >= 0){
            path1 <- rbind(path1,StCoordLine(rtp[1],rtp[2],wulff=wulff))
        } else {
            path2 <- rbind(path2,StCoordLine(rtp[1],rtp[2],wulff=wulff))
        }
    }
    list(path1=path1,path2=path2)
}

stereonet(trd=c(pi/3,pi/4),plg=c(pi/4,pi/5),option=1,
          wulff=FALSE,show.grid=TRUE,lwd=2,pch=19)
stereonet(trd=c(pi/3,pi/4),plg=c(pi/4,pi/5),option=2,add=TRUE,
          wulff=FALSE,show.grid=FALSE,lwd=2,pch=20)
stereonet(trd=pi/3,plg=pi/4,coneAngle=pi/10,option=4,add=TRUE,
          wulff=FALSE,show.grid=FALSE,lwd=2)
