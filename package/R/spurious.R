#' @title Spurious correlation
#' @description Calculate the 'null correlation' of ratios
#' @details Implements the spurious correlation formula of Pearson (1897)
#' @param mw the mean of variable \code{w}
#' @param mx the mean of variable \code{x}
#' @param my the mean of variable \code{y}
#' @param mz the mean of variable \code{z}
#' @param sw the standard deviation of variable \code{w}
#' @param sx the standard deviation of variable \code{x}
#' @param sy the standard deviation of variable \code{y}
#' @param sz the standard deviation of variable \code{z}
#' @param rwx the correlation coefficient between \code{w} and \code{x}
#' @param rwy the correlation coefficient between \code{w} and \code{y}
#' @param rwz the correlation coefficient between \code{w} and \code{z}
#' @param rxy the correlation coefficient between \code{x} and \code{y}
#' @param rxz the correlation coefficient between \code{x} and \code{z}
#' @param ryz the correlation coefficient between \code{y} and \code{z}
#' @return the null correlation coefficient
#' @rdname rwyxz
#' @examples
#' rxzyz(mx=100,my=100,mz=100,sx=1,sy=1,sz=10)
#' @export
rwyxz <- function(mw,mx,my,mz,sw,sx,sy,sz,rwx=0,rwy=0,rwz=0,rxy=0,rxz=0,ryz=0){
    out <- (
        rwx*(sw/mw)*(sx/mx) - rwz*(sw/mw)*(sz/mz) -
        rxy*(sx/mx)*(sy/my) + ryz*(sy/my)*(sz/mz)
    )/(
        sqrt((sw/mw)^2+(sy/my)^2-2*rwy*(sw/mw)*(sy/my)) *
        sqrt((sx/mx)^2+(sz/mz)^2-2*rxz*(sx/mx)*(sz/mz))
    )
    out
}
#' @rdname rwyxz
#' @export
ryxy <- function(mx,my,sx,sy,rxy=0){
    rwyxz(mw=my,mx=mx,my=1,mz=my,sw=sy,sx=sx,sy=0,sz=sy,
          rwx=rxy,rwy=0,rwz=1,rxy=0,rxz=rxy,ryz=0)
}
#' @rdname rwyxz
#' @export
rxzyz <- function(mx,my,mz,sx,sy,sz,rxy=0,rxz=0,ryz=0){
    rwyxz(mw=my,mx=mx,my=my,mz=mz,sw=sy,sx=sx,sy=sz,sz=sz,
          rwx=rxy,rwy=ryz,rwz=ryz,rxy=rxz,rxz=rxz,ryz=1)
}
