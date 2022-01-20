idir <- "~/Documents/Programming/R/geostats/quizzes/"
odir <- "~/Desktop/GEOL0061/"
rmarkdown::render(input=paste0(idir,"extra.Rmd"),
                  output_file=paste0(odir,"extra.pdf"))

if (TRUE){ # parse
    pages <- c(2,4,7,8,9,11,14,17,18,21,23,25,27)
    np <- length(pages)
    set.seed(1)
    if (TRUE){ # randomise
        NN <- round(runif(n=np,min=0,max=20))
    } else {
        NN <- rep("",np)
    }
    for (i in 1:(np-1)){
        cmd <- paste0("pdftk ",odir,"extra.pdf cat ",
                      pages[i],"-",pages[i+1]-1,
                      " output ",odir,"q",i,NN[i],".pdf")
        system(cmd)
    }
}
