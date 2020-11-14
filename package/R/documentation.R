#' pH data
#'
#' pH measurements in 20 samples of rain water
#' 
#' @name pH
#' @docType data
#' @keywords data
#' @examples
#' data(pH,package='geostats')
#' hist(pH)
NULL

#' clast size data
#'
#' 20 clast size measurements, in cm
#' 
#' @name clasts
#' @docType data
#' @keywords data
#' @examples
#' data(clasts,package='geostats')
#' d <- density(log(clasts))
#' plot(d)
NULL

#' porosity data
#'
#' 20 porosity measurements, as fractions
#' 
#' @name porosity
#' @docType data
#' @keywords data
#' @examples
#' data(porosity,package='geostats')
#' plot(density(logit(porosity)))
NULL

#' detrital zircon U-Pb data
#'
#' detrital zircon U-Pb data of 5 sand samples from China
#' 
#' @name DZ
#' @docType data
#' @keywords data
#' @examples
#' data(DZ,package='geostats')
#' qqplot(DZ[['Y']],DZ[['5']])
NULL

#' Rb-Sr data
#'
#' synthetic dataset of 8 Rb-Sr analysis that form a 1Ga isochron
#' 
#' @name rbsr
#' @docType data
#' @keywords data
#' @examples
#' data(rbsr,package='geostats')
#' plot(rbsr[,'RbSr'],rbsr[,'SrSr'])
#' fit <- lm(SrSr ~ RbSr,data=rbsr)
#' abline(fit)
NULL

#' declustered earthquake data
#'
#' dataset of 28267 earthquakes between 1769 and 2016, with
#' aftershocks and precursor events removed
#' 
#' @name declustered
#' @docType data
#' @keywords data
#' @references Mueller, C.S., 2019. Earthquake catalogs for the USGS
#'     national seismic hazard maps. Seismological Research Letters,
#'     90(1), pp.251-261.
#' @examples
#' data(declustered,package='geostats')
#' quakesperyear <- countQuakes(declustered,minmag=5.0,from=1917,to=2016)
#' table(quakesperyear)
NULL

#' earthquake data
#'
#' dataset of 20000 earthquakes between 2017 and 2000, downloaded from
#' the USGS earthquake database
#' (\url{https://earthquake.usgs.gov/earthquakes/search/}).
#' 
#' @name earthquakes
#' @docType data
#' @keywords data
#' @examples
#' data(earthquakes,package='geostats')
#' gutenberg(earthquakes$mag)
NULL

#' Finnish lake data
#'
#' Table of 2327 Finnish lakes, extracted from a hydroLAKES database.
#' 
#' @name Finland
#' @docType data
#' @keywords data
#' @references Lehner, B., and Doll, P. (2004), Development and
#'     validation of a global database of lakes, reservoirs and
#'     wetlands, Journal of Hydrology, 296(1), 1-22, doi:
#'     10.1016/j.jhydrol.2004.03.028.
#' @examples
#' data(Finland,package='geostats')
#' sf <- sizefrequency(Finland$area)
#' size <- sf[,'size']
#' freq <- sf[,'frequency']
#' plot(size,freq,log='xy')
#' fit <- lm(log(freq) ~ log(size))
#' lines(exp(predict(fit)))
NULL

#' foram count data
#'
#' Planktic foraminifera counts in surface sediments in the Atlantic ocean.
#' 
#' @name forams
#' @docType data
#' @keywords data
#' @examples
#' data(forams,package='geostats')
#' abundant <- forams[,c('quinqueloba','pachyderma','incompta',
#'                       'glutinata','bulloides')]
#' other <- rowSums(forams[,c('uvula','scitula')])
#' dat <- cbind(abundant,other)
#' chisq.test(dat)
NULL

#' world population
#'
#' The world population from 1750 until 2014
#' 
#' @name pop
#' @docType data
#' @keywords data
#' @examples
#' data(pop,package='geostats')
#' plot(pop)
NULL

#' British coast
#'
#' a 512 x 512 pixel image of the British coast line
#' 
#' @name Britain
#' @docType data
#' @keywords data
#' @examples
#' data(Britain,package='geostats')
#' p <- par(mfrow=c(1,2))
#' image(Britain)
#' fractaldim(Britain)
NULL

#' rivers on Corsica
#'
#' a 512 x 512 pixel image of the river network on Corsica
#' 
#' @name Corsica
#' @docType data
#' @keywords data
#' @examples
#' data(Corsica,package='geostats')
#' p <- par(mfrow=c(1,2))
#' image(Corsica)
#' fractaldim(Corsica)
NULL

#' fractures
#'
#' a 512 x 512 pixel image of a fracture network
#' 
#' @name fractures
#' @docType data
#' @keywords data
#' @examples
#' data(fractures,package='geostats')
#' p <- par(mfrow=c(1,2))
#' image(fractures)
#' fractaldim(fractures)
NULL

#' A-CN-K compositions
#'
#' Synthetic A (Al2O3), CN (CaO+Na2O), K (K2O) data table
#' 
#' @name ACNK
#' @docType data
#' @keywords data
#' @examples
#' data(ACNK,package='geostats')
#' ternary(ACNK,type='p',labels=c(expression('Al'[2]*'O'[3]),
#'                                expression('CaO+Na'[2]*'O'),
#'                                expression('K'[2]*'O')))
NULL

#' composition of Namib dune sand
#'
#' major element compositions of 16 Namib sand samples
#' 
#' @name major
#' @docType data
#' @keywords data
#' @examples
#' data(major,package='geostats')
#' comp <- clr(major)
#' pc <- prcomp(comp)
#' biplot(pc)
NULL

#' composition of oceanic basalts
#'
#' major element compositions of 227 island arc basalts (IAB), 221 mid
#' oceanic ridge basalts (MORB) and 198 ocean island basalts
#' (OIB). This dataset can be used to train supervised learning
#' algorithms.
#' 
#' @name training
#' @docType data
#' @keywords data
#' @examples
#' library(MASS)
#' data(training,package='geostats')
#' qd <- qda(affinity ~ ., data=training)
#' pr <- predict(qd)
#' table(training$affinity,pr$class)
NULL

#' composition of oceanic basalts
#'
#' major element compositions of 64 island arc basalts (IAB), 23 mid
#' oceanic ridge basalts (MORB) and 60 ocean island basalts
#' (OIB). This dataset can be used to test supervised learning
#' algorithms.
#' 
#' @name test
#' @docType data
#' @keywords data
#' @examples
#' library(MASS)
#' data(training,package='geostats')
#' data(test,package='geostats')
#' qd <- qda(affinity ~ ., data=training)
#' pr <- predict(qd,newdata=test[,-1])
#' table(test$affinity,pr$class)
NULL

#' A-F-M data
#'
#' (Na2O + K2O) - FeO - MgO compositions of 630 calc-alkali basalts
#' from the Cascade Mountains and 474 tholeiitic basalts from Iceland.
#' 
#' @name AFM
#' @docType data
#' @keywords data
#' @examples
#' data(AFM,package='geostats')
#' ternary(AFM[,-1])
NULL
