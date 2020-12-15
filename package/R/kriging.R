getsnr <- function(lsnr){
    s <- exp(lsnr[1]) # sill
    n <- s*exp(lsnr[2])/(1+exp(lsnr[2])) # nugget
    r <- exp(lsnr[3]) # range
    c(s,n,r)
}

semivarmod <- function(h,lsnr,model='spherical'){
    snr <- getsnr(lsnr)
    s <- snr[1]
    n <- snr[2]
    r <- snr[3]
    out <- rep(n,length(h))
    i <- (h>0 & h<r)
    if (model=='spherical'){
        out[h==0] <- n
        out[i] <- n + (s-n)*(1.5*h[i]/r-0.5*(h[i]/r)^3)
        out[h>=r] <- s
    } else {
        stop('Invalid semivariogram model')
    }
    out
}

semivariogram <- function(x,y,z,bw=NULL,plot=TRUE,fit=FALSE){
    d <- as.matrix(dist(cbind(x,y)))
    if (is.null(bw)){
        N <- length(z)
        q4 <- min(1,4/N)
        bw <- quantile(d,probs=q4)
        nq <- 10
        sv <- rep(0,nq)
        for (i in 1:nq){
            ij <- which((d>=((i-1)*bw)) & d< (i*bw),arr.ind=TRUE)
            sv[i] <- mean((z[ij[,1]]-z[ij[,2]])^2)/2
        }
    }
    h <- bw*((1:nq)-1/2) # lag
    if (plot){
        plot(h,sv,type='p',xlab='lag',ylab='var',ylim=c(0,max(sv)))
    }
    if (fit){
        misfit <- function(lsnr,h,sv){
            svp <- semivarmod(h,lsnr)
            sum((sv-svp)^2)
        }
        lsnr <- optim(par=c(log(max(sv)),-2,log(max(h))),
                      fn=misfit,h=h,sv=sv)$par
        if (plot){
            x <- seq(from=0,to=max(h),length.out=20)
            y <- semivarmod(h=x,lsnr=lsnr)
            lines(x,y)
        }
        out <- getsnr(lsnr)
        names(out) <- c('sill','nugget','range')
    } else {
        out <- list(h=h,sv=sv)
    }
    invisible(out)
}


data('meuse',package='sp')
dat <- meuse
fit <- semivariogram(dat$x,dat$y,log(dat$cadmium),fit=TRUE)

