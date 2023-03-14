library(geostats)

X <- meuse$x # Easting
Y <- meuse$y # Northing
Z <- log(meuse$zinc)
colourplot(X=X,Y=Y,Z=Z,key.title=title('ln[Zn]'))

svm <- semivariogram(x=X,y=Y,z=Z)

svm_exp <- semivariogram(x=X,y=Y,z=Z,model='exponential')

xi <- seq(from=min(X),to=max(X),length.out=50)
yi <- seq(from=min(Y),to=max(Y),length.out=50)
zi <- kriging(x=X,y=Y,z=Z,svm=svm,xi=xi,yi=yi,grid=TRUE)
colourplot(x=xi,y=yi,z=zi,key.title=title('ln[Zn]'))

zi <- kriging(x=X,y=Y,z=Z,svm=svm,xi=xi,yi=yi,grid=TRUE,err=TRUE)
colourplot(x=xi,y=yi,z=sqrt(zi),key.title=title('s[Zn]/Zn'))