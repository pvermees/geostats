# the 'print' function is needed to show intermediate
# results when running commands from a .R file
print(pi)

toss <- function(){
  r <- runif(1) # create a random number between 0 and 1
  if (r<0.5){
  print("head")
  } else {
    print("tail")
  }
}

fibonnaci <- function(n=5){ # 5 is the default value
  if (n < 3) { stop('n must be at least 3') }
  # seed the output vector with 0 and 1:
  s <- c(0,1)
  # loop through all numbers from 3 to n:
  for (i in 3:n){
    s[i] <- s[i-1] + s[i-2]
  }
  return(s)
}

library(geostats)

s <- sierpinski()
image(s,axes=FALSE)
