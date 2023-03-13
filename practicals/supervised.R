library(geostats)
library(MASS)

ld <- lda(Species ~ ., data=iris)
newflower <- data.frame(Sepal.Length=6.0,Sepal.Width=3.0,
Petal.Length=5.0,Petal.Width=1.5)
pred <- predict(ld,newdata=newflower)

qd <- qda(Species ~ ., data=iris)

library(rpart)
tree <- rpart(Species ~ ., data=iris, method="class")
plot(tree)
text(tree)

text(tree, use.n=TRUE)