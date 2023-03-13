library(geostats)

X <- rbind(c(-1,7),c(3,2),c(4,3))
colnames(X) <- c('a','b')
PCA2D(X)

d <- dist(X)           # create a Euclidean distance matrix
conf <- cmdscale(d)    # classical MDS
plot(conf,type='n')    # create an empty plot
text(conf,labels=1:3)  # add text labels to the empty plot

pc <- prcomp(USArrests,scale.=TRUE)
biplot(pc)

conf <- cmdscale(eurodist)
plot(conf,type='n',asp=1)
text(conf,labels=labels(eurodist))

library(MASS)
mds <- isoMDS(eurodist)
conf <- mds$points
plot(conf,type='n',asp=1)
text(conf,labels=labels(eurodist))

ylim <- rev(range(conf[,2])) # reverse the minimum and maximum values
plot(conf,type='n',asp=1,ylim=ylim) # change the y-axis limits

sh <- Shepard(d=eurodist,x=conf)
stress <- signif(mds$stress,2)
plot(sh,main=paste0('stress=',stress))

measurements <- iris[,-5]
species <- iris[,5]
plot(measurements,pch=as.numeric(species))

fit <- kmeans(measurements,centers=3)

tree <- hclust(dist(measurements))
plot(tree)

treecut <- cutree(tree,k=3)
table(treecut,species)