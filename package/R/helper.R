mysign <- function(x){
    out <- rep(1,length(x))
    out[x<0] <- -1
    out
}
