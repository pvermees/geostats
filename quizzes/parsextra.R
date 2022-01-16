idir <- "~/Documents/Programming/R/geostats/quizzes/"
odir <- "~/Desktop/extra/"
rmarkdown::render(paste0(idir,"extra.Rmd"))
pages <- c(2,4,7,8,9,11,13,16,19,22)
np <- length(pages)
set.seed(1)
NN <- round(runif(n=nq,min=0,max=20))
for (i in 1:(np-1)){
    cmd <- paste0("pdftk ",idir,"extra.pdf cat ",
                  pages[i],"-",pages[i+1]-1,
                  " output ",odir,"q",i,NN[i],".pdf")
    system(cmd)
}
