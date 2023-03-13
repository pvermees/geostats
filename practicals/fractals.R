library(geostats)

sf <- sizefrequency(Finland$area)
plot(frequency~size,data=sf,log='xy')
fit <- lm(log(frequency)~log(size),data=sf)
lines(x=sf$size,y=exp(predict(fit)))

gutenberg(earthquakes$mag)

k <- koch(n=5)
fit <- fractaldim(k)