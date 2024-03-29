% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformations.R
\name{alr}
\alias{alr}
\title{additive logratio transformation}
\usage{
alr(dat, inverse = FALSE)
}
\arguments{
\item{dat}{an \eqn{n} column data frame or matrix}

\item{inverse}{if \code{TRUE}, applies the inverse alr
tranformation}
}
\value{
If \code{inverse=FALSE}, returns an \eqn{(n-1) \times m}
    matrix of logratios; otherwise returns an \eqn{(n+1) \times m}
    matrix of compositional data whose columns add up to 1.
}
\description{
Maps compositional data from an \eqn{n}-dimensional
    simplex to an (\eqn{n-1})-dimensional Euclidean space with
    Aitchison's additive logratio transformation.
}
\examples{
xyz <- rbind(c(0.03,99.88,0.09),
             c(70.54,25.95,3.51),
             c(72.14,26.54,1.32))
colnames(xyz) <- c('a','b','c')
rownames(xyz) <- 1:3
uv <- alr(xyz)
XYZ <- alr(uv,inverse=TRUE)
xyz/XYZ
}
