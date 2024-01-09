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

sink(file="tables.tex")
"\\chapter{Tables}\\label{ch:tables}"
xtable::xtable(normtab,caption='Upper quantiles of the standard normal distribution (lower quantiles are negative equivalents for reasons of symmetry)')
xtable::xtable(ttab,caption='Upper quantiles of the t-distribution with $df$ degrees of freedom (lower quantiles are negative equivalents for reasons of symmetry)')
xtable::xtable(x2tab,caption='Quantiles of the Chi-square distribution with $df$ degrees of freedom')
sink()
