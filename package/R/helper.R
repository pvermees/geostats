mysign <- function(x){
    out <- rep(1,length(x))
    out[x<0] <- -1
    out
}

# from eponymous fuction in pracma package
inpolygon <- function (x, y, xp, yp, boundary = TRUE){
    stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y), 
        is.numeric(xp), is.numeric(yp), length(xp) == length(yp))
    n <- length(x)
    np <- length(xp)
    if (xp[1] != xp[np] || yp[1] != yp[np]) {
        xp <- c(xp, xp[1])
        yp <- c(yp, yp[1])
        np <- np + 1
    }
    inpoly <- rep(FALSE, n)
    onpoly <- rep(FALSE, n)
    j <- np
    for (i in 1:np) {
        dxp <- xp[j] - xp[i]
        dyp <- yp[j] - yp[i]
        dist <- dxp * (y - yp[i]) - (x - xp[i]) * dyp
        idx1 <- (((yp[i] <= y & y < yp[j]) | (yp[j] <= y & y < 
            yp[i])) & (0 < dist * dyp))
        inpoly[idx1] <- !inpoly[idx1]
        idx2 <- (((yp[i] <= y & y <= yp[j]) | (yp[j] <= y & y <= 
            yp[i])) & ((xp[i] <= x & x <= xp[j]) | (xp[j] <= 
            x & x <= xp[i])) & (0 == dist | !dxp))
        onpoly[idx2] <- TRUE
        j <- i
    }
    if (boundary) inpoly[onpoly] <- TRUE
    else inpoly[onpoly] <- FALSE
    return(inpoly)
}

inhull <- function(x,y,xi,yi,buffer=0.05){
    dx <- buffer*(max(x)-min(x))
    dy <- buffer*(max(y)-min(y))
    X <- c(x-dx,x-dx,x+dx,x+dx)
    Y <- c(y-dx,y+dy,y-dy,y+dy)
    ch <- grDevices::chull(x=X,y=Y)
    xh <- X[ch]
    yh <- Y[ch]
    inpolygon(x=xi,y=yi,xp=xh,yp=yh)
}


#' @title colour plot
#' @description Adds a colour bar to a scatter plot and/or filled
#'     contour plot.  This function, which is based on base \code{R}'s
#'     \code{filled.contour} function, is useful for visualising
#'     kriging results.
#' @param x numerical vector of \eqn{n} equally spaced values to be
#'     used in the contour plot.
#' @param y numerical vector of \eqn{m} equally spaced values to be
#'     used in the contour plot.
#' @param z an \eqn{n \times m} matrix of numerical values to be used
#'     in the contour plot.
#' @param X numerical vector of \eqn{N} values to be used in the
#'     scatter plot.
#' @param Y numerical vector of \eqn{N} values to be used in the
#'     scatter plot.
#' @param Z numerical vector of \eqn{N} values to be used in the
#'     scatter plot.
#' @param nlevels number of levels to be used in the contour plot.
#' @param colspec colour specification (e.g., \code{rainbow},
#'     \code{grey.colors}, \code{heat.colors}, \code{topo.colors}).
#' @param pch plot character (\code{21} -- \code{25}).
#' @param cex plot character magnification.
#' @param plot.title statements that add titles to the main plot.
#' @param plot.axes statements that draw axes on the main plot. This
#'     overrides the default axes.
#' @param key.title statements that add titles for the plot key.
#' @param key.axes statements that draw axes on the plot key.  This
#'     overrides the default axis.
#' @param asp the y/x aspect ratio, see \code{plot.window}.
#' @param xaxs the x axis style.  The default is to use internal
#'     labelling.
#' @param yaxs the y axis style.  The default is to use internal
#'     labelling.
#' @param las the style of labelling to be used.  The default is to
#'     use horizontal labelling.
#' @param axes logicals indicating if axes should be drawn.
#' @param frame.plot logicals indicating if a box should be drawn, as
#'     in \code{plot.default}.
#' @param extra (optional) extra intructions to be carried out in the
#'     main plot window, such as text annotations.
#' @param ... additional graphical parameters
#' @return no return value
#' @import grDevices
#' @import graphics
#' @examples
#' data('meuse',package='geostats')
#' colourplot(X=meuse$x,Y=meuse$y,Z=log(meuse$zinc))
#' @export
colourplot <- function (x, y, z, X, Y, Z, nlevels=20, colspec=hcl.colors,
                        pch = 21, cex = 1, plot.title, plot.axes, key.title,
                        key.axes, asp = NA, xaxs = "i", yaxs = "i",
                        las = 1, axes = TRUE, frame.plot = axes, extra, ...) {
    cnts <- !(missing(x)|missing(y)|missing(z))
    pnts <- !(missing(X)|missing(Y)|missing(Z))
    if (cnts){
        xX <- x
        yY <- y
        zZ <- z
        if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
            stop("increasing 'x' and 'y' values expected")
    } else {
        xX <- NULL
        yY <- NULL
        zZ <- NULL     
    }
    if (pnts){
        xX <- c(xX,X)
        yY <- c(yY,Y)
        zZ <- c(zZ,Z)
    }
    if (!pnts & !cnts){
        stop('Both x,y,z and X,Y,Z are missing.')
    }
    xlim <- range(xX,finite=TRUE)
    ylim <- range(yY,finite=TRUE)
    zlim <- range(zZ,finite=TRUE)
    # add a bit of padding
    dx <- diff(xlim)
    dy <- diff(ylim)
    dz <- diff(zlim)
    xlim <- xlim + dx*c(-1,1)/20
    ylim <- ylim + dy*c(-1,1)/20
    zlim <- zlim + dz*c(-1,1)/100
    levels <- seq(from=zlim[1],to=zlim[2],length.out=nlevels)
    dl <- diff(zlim)
    levcol <- colspec(nlevels)
    if (pnts){
        ci <- ceiling(nlevels*(Z-levels[1])/dl)
        ptscol <- levcol[ci]
    }
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = levcol)
    if (missing(key.axes)) {
        if (axes) 
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (cnts){
        .filled.contour(x, y, z, levels, levcol)
    }
    if (pnts){
        points(X, Y, pch=pch, cex=cex, bg=ptscol)
    }
    if (!missing(extra)){
        extra
    }
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            Axis(xlim, side = 1)
            Axis(ylim, side = 2)
        }
    }
    else plot.axes
    if (frame.plot) 
        box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
}
