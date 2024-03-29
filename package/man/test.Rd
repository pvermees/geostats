% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentation.R
\docType{data}
\name{test}
\alias{test}
\title{composition of a further 147 oceanic basalts}
\description{
Major element compositions of 64 island arc basalts (IAB), 23 mid
oceanic ridge basalts (MORB) and 60 ocean island basalts
(OIB). This dataset can be used to test supervised learning
algorithms.
}
\examples{
library(MASS)
data(training,package='geostats')
ld <- lda(x=alr(training[,-1]),grouping=training[,1])
data(test,package='geostats')
pr <- predict(ld,newdata=alr(test[,-1]))
table(test$affinity,pr$class)
}
\references{
Vermeesch, P. ``Tectonic discrimination diagrams
    revisited.''  Geochemistry, Geophysics, Geosystems 7.6 (2006).
}
\keyword{data}
