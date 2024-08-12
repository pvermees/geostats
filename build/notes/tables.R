p <- c(0.55,0.6,0.65,0.7,0.75,0.8,signif(pnorm(1),3),
       0.9,0.95,0.975,signif(pnorm(2),3),0.99)
np <- length(p)
normtab <- matrix(qnorm(p),nrow=1)
colnames(normtab) <- p

nr <- 20
dfs <- 1:nr
ttab <- matrix(NA,nrow=nr,ncol=np)
colnames(ttab) <- p
rownames(ttab) <- dfs
for (df in dfs){
    ttab[df,] <- qt(p,df)
}

p <- c(0.025,0.05,seq(from=0.1,to=0.9,by=0.1),0.95,0.975)
np <- length(p)
x2tab <- matrix(NA,nrow=nr,ncol=np)
colnames(x2tab) <- p
rownames(x2tab) <- dfs
for (df in dfs){
    x2tab[df,] <- qchisq(p,df)
}

wtab <- function(alpha){
    n2 <- 3:20
    out <- matrix(NA,nrow=length(n2),ncol=max(n2))
    rownames(out) <- n2
    colnames(out) <- 1:max(n2)
    for (i in seq_along(n2)){
        for (n1 in 1:n2[i]){
            scaled <- qwilcox(alpha,n1,n2[i],lower.tail=TRUE)
            if (scaled>0){ # unscale
                out[i,n1] <- round(scaled + n1*(n1+1)/2 - 1)
            }
        }
    }
    out[out<0] <- NA
    out
}    

sink(file="tables.tex")
"\\chapter{Tables}\\label{ch:tables}"
xtable::xtable(normtab,caption='Upper quantiles of the standard normal distribution (lower quantiles are negative equivalents for reasons of symmetry)')
xtable::xtable(ttab,caption='Upper quantiles of the t-distribution with $df$ degrees of freedom (lower quantiles are negative equivalents for reasons of symmetry)')
xtable::xtable(x2tab,caption='Quantiles of the Chi-square distribution with $df$ degrees of freedom')
xtable::xtable(wtab(alpha=0.025),caption='Critical values (upper bounds of the rejection region) for a Wilcoxon test for $\\alpha=0.025$. Column labels ($n_1$) and row labels ($n_2$) are the sample sizes, with ${n_1}\\leq{n_2}$.',digits=0)
xtable::xtable(wtab(alpha=0.05),caption='Critical values (upper bounds of the rejection region) for a Wilcoxon test for $\\alpha=0.05$. Column labels ($n_1$) and row labels ($n_2$) are the sample sizes, with ${n_1}\\leq{n_2}$.',digits=0)
sink()
