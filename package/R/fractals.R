#' @title Sierpinski carpet
#' @description returns a matrix of 0s and 1s that form a Sierpinski
#'     fractal.
#' @details The Sierpinski carpet is two dimensional fractal, which is
#'     generated using a recursive algorithm that is built on a grid
#'     of eight black squares surrounding a white square. Each level
#'     of recursion replaces each black square by the same pattern.
#' @param n an integer value controling the number of recursive
#'     levels.
#' @return a square matrix with 0s and 1s.
#' @examples
#' g <- sierpinski(n=5)
#' image(g,col=c('white','black'),axes=FALSE,asp=1)
#' @export
sierpinski <- function(n=5){
    X <- matrix(1,1,1)
    for (i in 1:n){
        Y <- cbind(X,X,X)
        Z <- cbind(X,0*X,X)
        X <- rbind(Y,Z,Y)
    }
    return(X)
}

#' @title Cantor set
#' @description Calculates or plots a Cantor set of fractal lines.
#' @details The Cantor set is generated using a recursive algorithm
#'     that is built on a line segment whose middle third is
#'     removed. Each level of recursion replaces each black line by
#'     the same pattern.
#' @param n an integer value controling the number of recursive
#'     levels.
#' @param plot logical.  If \code{TRUE}, the Cantor set is plotted,
#'     otherwise a list of breaks and counts is returned.
#' @param add logical (only used if \code{plot=TRUE}). If
#'     \code{add=FALSE}, then a brand new figure is created; otherwise
#'     the Cantor set is added to an existing plot.
#' @param Y y-value for the plot (only used if \code{plot=TRUE}).
#' @param lty line type (see \code{pars()} for details)
#' @param col colour of the Cantor lines.
#' @param ... optional arguments to be passed on to \code{matplot} or
#'     \code{matlines}.
#' @return a square matrix with 0s and 1s.
#' @examples
#' g <- sierpinski(n=5)
#' image(g,col=c('white','black'),axes=FALSE,asp=1)
#' @export
cantor <- function(n=5,plot=FALSE,add=FALSE,Y=0,
                   lty=1,col='black',...){
    if (n <= 0){
        x <- c(0,1)
    } else {
        X <- cantor(n-1,Y=Y)
        x <- (1/3)*c(X,X+2)
    }
    if (plot){
        xx <- rbind(x[1:(length(x)-1)],x[-1])
        yy <- matrix(Y,nrow(xx),ncol(xx))
        i <- seq(from=1,to=ncol(xx),by=2)
        if (add){
            matlines(xx[,i],yy[,i],lty=lty,col=col,...)
        } else {
            matplot(xx[,i],yy[,i],type='l',lty=lty,col=col,...)
        }
    }
    invisible(x)
}

#' @title Koch snowflake
#' @description Calculates or plots a Koch set of fractal lines.
#' @details The Koch set is generated using a recursive algorithm that
#'     is built on a triangular hat shaped line segment. Each level of
#'     recursion replaces each linear segment by the same pattern.
#' @param n an integer value controling the number of recursive
#'     levels.
#' @param plot logical.  If \code{TRUE}, the Koch flake is plotted.
#' @param res the number of pixels in each side of the output matrix
#' @return a \code{res x res} matrix with 0s and 1s
#' @examples
#' k <- koch(n=5)
#' d <- fractaldim(k,plot=FALSE)
#' print(d)
#' @export
koch <- function(n=4,plot=TRUE,res=512){
    i <- function(x,res=512){
        max(1,round(x*res/100))
    }
    recurse <- function(n=n,x1=0,y1=100,x5=0,y5=0,m=matrix(0,res,res)){
        if (n == 0){
            m[i(y1,res=res),i(x1,res=res)] <- 1
            m[i(y5,res=res),i(x5,res=res)] <- 1
            if (plot) lines(x=c(y1,y5),y=c(x1,x5))
        } else {
            deltaX = x5 - x1
            deltaY = y5 - y1
            x2 = x1 + deltaX / 3
            y2 = y1 + deltaY / 3
            x3 = (0.5 * (x1+x5) + sqrt(3) * (y1-y5)/6)
            y3 = (0.5 * (y1+y5) + sqrt(3) * (x5-x1)/6)
            x4 = x1 + 2 * deltaX /3
            y4 = y1 + 2 * deltaY /3
            m <- recurse(n=n-1,x1=x1,y1=y1,x5=x2,y5=y2,m=m)
            m <- recurse(n=n-1,x1=x2,y1=y2,x5=x3,y5=y3,m=m)
            m <- recurse(n=n-1,x1=x3,y1=y3,x5=x4,y5=y4,m=m)
            m <- recurse(n=n-1,x1=x4,y1=y4,x5=x5,y5=y5,m=m)
        }
        m
    }
    if (plot) plot(c(0,100),c(0,100),type='n',axes=FALSE,xlab='',ylab='')
    height <- 100*3/4
    width <- 100
    xStart <- width/2 - height/2
    dx <- res/100
    m <- recurse(n=n,x1=xStart+dx,y1=height-dx,x5=xStart+height-dx,y5=height-dx)
    m <- recurse(n=n,x1=xStart+height-dx,y1=height-dx,x5=xStart+height/2,y5=dx,m=m)
    m <- recurse(n=n,x1=xStart+height/2,y1=dx,x5=xStart+dx,y5=height-dx,m=m)
    m
}

#' @title box counting
#' @description count the number of boxes needed to cover all the 1s
#'     in a matrix of 0s and 1s.
#' @param mat a square square matrix of 0s and 1s. Must be a power of 2.
#' @param size the size (pixels per side) of the boxes. Should be a power of 2.
#' @param plot logical. If \code{TRUE}, plots the results on a log-log
#'     scale.
#' @examples
#' g <- sierpinski(n=5)
#' boxcount(mat=g,size=16)
#' @export
boxcount <- function(mat,size){
    n2 <- floor(log2(min(dim(mat))))
    n <- 2^(n2-log2(size)-1)
    nboxes <- n^2
    for(j in 1:n){
        row_id <- ((j-1)*size+1):(j*size)
        for(k in 1:n){
            col_id <- ((k-1)*size+1):(k*size)
            if(sum(mat[row_id,col_id]) == 0){
                nboxes <- nboxes - 1
            }
        }
    }
    nboxes
}

#' @title calculate the fractal dimension
#' @description performs box counting on a matrix of 0s and 1s.
#' @param mat a square matrix of 0s and 1s. Size must be a power of 2.
#' @param plot logical. If \code{TRUE}, plots the results on a log-log
#'     scale.
#' @param ... optional arguments to the generic \code{points} function.
#' @examples
#' g <- sierpinski(n=5)
#' fractaldim(g)
#' @export
fractaldim <- function(mat,plot=TRUE,...){
    n2 <- floor(log2(min(dim(mat))))
    n <- rep(2,n2)^(0:(n2-1))
    size <- rev(n)
    nboxes <- n^2
    for(i in 1:length(size)){
        nboxes[i] <- boxcount(mat=mat,size=size[i])
    }
    fit <- lm(log(nboxes) ~ log(size))
    if (plot){
        plot(log(nboxes) ~ log(size),type='n',
             bty='n',xlab='ln[size of boxes]',
             ylab=expression('ln[number of boxes]'))
        abline(fit)
        points(log(nboxes) ~ log(size),...)
        legend('topright',bty='n',
               legend=paste0('y = ',signif(fit$coef[1],3),
                             signif(fit$coef[2],3),' x'))
    }
    invisible(fit)
}

#' @title count the number of earthquakes per year
#' @description counts the number of earthquakes per year that fall
#'     between two magnitude limits
#' @param qdat a data frame containing columns named \code{mag} and
#'     \code{year}.
#' @param minmag minimum magnitude
#' @param from first year
#' @param to last year
#' @examples
#' data(declustered,package='geostats')
#' quakesperyear <- countQuakes(declustered,minmag=5.0,from=1917,to=2016)
#' table(quakesperyear)
#' @export
countQuakes <- function(qdat,minmag,from,to){
    bigenough <- (declustered$mag >= minmag)
    youngenough <- (declustered$year >= from)
    oldenough <- (declustered$year <= to)
    goodenough <- (bigenough & youngenough & oldenough)
    table(qdat$year[goodenough])
}

#' @title create a Gutenberg-Richter plot
#' @description calculate a semi-log plot with earthquake magnitude on
#'     the horizontal axis,and the cumulative number of earthquakes
#'     exceeding any given magnitude on the vertical axis.
#' @param mag a vector of earthquake magnitudes
#' @param n the number of magnitudes to evaluate
#' @param ... optional arguments to the generic \code{points} function.
#' @examples
#' data(declustered,package='geostats')
#' gutenberg(declustered$mag)
#' @export
gutenberg <- function(mag,n=10,...){
    sf <- sizefrequency(mag,n=n,log=FALSE)
    X <- sf[,'size']
    Y <- log10(sf[,'frequency']/length(mag))
    plot(x=X,y=Y,type='n',bty='n',xlab='magnitude',
         ylab=expression('log'[10]*'[N/N'[o]*']'))
    fit <- lm(Y ~ X)
    lines(X,fit$coef[1]+fit$coef[2]*X)
    points(x=X,y=Y,...)
    legend('topright',paste0('y = ',signif(fit$coef[1],3),
                             signif(fit$coef[2],3),' x'),bty='n')
}

#' @title calculate the size-frequency distribution of things
#' @description calculate the number of items exceeding a certain size
#' @param dat a numerical vector
#' @param n the number of sizes to evaluate
#' @param log logical. If \code{TRUE}, uses a log spacing for the
#'     sizes at which the frequencies are evaluated
#' @examples
#' data(Finland,package='geostats')
#' sf <- sizefrequency(Finland$area)
#' size <- sf[,'size']
#' freq <- sf[,'frequency']
#' plot(size,freq,log='xy')
#' fit <- lm(log(freq) ~ log(size))
#' lines(x=size,y=exp(predict(fit)))
#' @export
sizefrequency <- function(dat,n=10,log=TRUE){
    if (log) d <- log(dat)
    else d <- dat
    m <- quantile(d,0.01)
    M <- quantile(d,0.99)
    size <- seq(from=m,to=M,length.out=n)
    freq <- rep(0,n)
    for (i in 1:n){
        freq[i] <- sum(d>size[i])
    }
    if (log) size <- exp(size)
    out <- cbind(size,freq)
    colnames(out) <- c('size','frequency')
    out
}
