setwd('~/Desktop/geostats/')
library(geostats)

pH <- catchments$pH

mean(pH)

sd(pH)

median(pH)

mad(pH,constant=1)

IQR(pH,type=1)

agetab <- table(catchments$age)
mod <- which.max(agetab)

dens <- density(pH)
mod <- dens$x[which.max(dens$y)]

skew <- function(x){
  mean((x-mean(x))^3)/sd(x)^3
}

skew(catchments$CaMg)

boxplot(catchments$CaMg)
