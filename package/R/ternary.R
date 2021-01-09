#' @title get x,y plot coordinates of ternary data
#' @description Helper function to generate bivariate plot coordinates
#'     for ternary data.
#' @param xyz an \code{n x 3} matrix or data frame
#' @return an \code{n x 2} numerical matrix
#' @examples
#' xyz <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,0,0))
#' xy <- xyz2xy(xyz)
#' plot(xy,type='l',bty='n')
#' @export
xyz2xy <- function(xyz){
    if (hasClass(xyz,'matrix','data.frame')){
        n <- nrow(xyz)
        x <- xyz[,1]
        y <- xyz[,2]
        z <- xyz[,3]
    } else {
        n <- 1
        x <- xyz[1]
        y <- xyz[2]
        z <- xyz[3]
    }
    xy <- matrix(0,nrow=n,ncol=2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- sin(pi/3)*x/(x+y+z)
    return(xy)
}

#' @title ternary diagrams
#' @description Plot points, lines or text on a ternary diagram.
#' @param xyz an \code{n x 3} matrix or data frame
#' @param f a three-element vector of multipliers for \code{xyz}
#' @param labels the text labels for the corners of the ternary
#'     diagram
#' @param add if \code{TRUE}, adds information to an existing ternary
#'     diagram
#' @param type one of \code{'n'} (empty plot), \code{'p'} (points),
#'     \code{'l'} (lines) or \code{'t'} (text).
#' @param ... optional arguments to the \code{points}, \code{lines} or
#'     \code{text} functions.
#' @examples
#' data(ACNK,package='geostats')
#' ternary(ACNK,type='p',labels=c(expression('Al'[2]*'O'[3]),
#'                                expression('CaO+Na'[2]*'O'),
#'                                expression('K'[2]*'O')))
#' @export
ternary <- function(xyz=NULL,f=rep(1,3),labels,
                    add=FALSE,type='p',...){
    if (!hasClass(xyz,'matrix','data.frame')){
        xyz <- matrix(xyz,nrow=1)
    }
    if (!add){
        if (missing(labels)){
            cnames <- colnames(xyz)
            if (is.null(cnames)){
                labels <- c('X','Y','Z')
            } else {
                labels <- cnames[1:3]
            }
        }
        corners <- rbind(c(1,0,0),c(0,1,0),c(0,0,1),c(1,0,0))
        xy <- xyz2xy(corners)
        graphics::plot(xy,type='l',asp=1,axes=FALSE,
                       ann=FALSE,bty='n')
        position <- c(3,1,1)
        for (i in 1:3){
            if (f[i]==1) lab <- labels[i]
            else lab <- paste0(f[i],'x',labels[i])
            graphics::text(xy[i,,drop=FALSE],labels=lab,pos=position[i],xpd=NA)
        }
    }
    XYZ <- t(apply(as.matrix(xyz,ncol=3),MARGIN=1,FUN='*',f))
    XYZ <- apply(XYZ,MARGIN=2,FUN='/',rowSums(XYZ))
    if (!is.null(xyz)){
        if (type=='p'){
            graphics::points(xyz2xy(XYZ),...)
        } else if (type=='l'){
            graphics::lines(xyz2xy(XYZ),...)
        } else if (type=='t'){
            graphics::text(xyz2xy(XYZ),labels=rownames(xyz),...)
        }
    }
}
