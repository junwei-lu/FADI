args=(commandArgs(TRUE))
##args=c(cohort_i, omega_l, p)
## m = 100, L = 80, p = 50 
i <- as.numeric(args[1])
l <- as.numeric(args[2])
p <- as.numeric(args[3])
.libPaths('/users/shutingshen/apps/R_inter')
library(gdsfmt)
gfile <- openfn.gds("/users/shutingshen/fast_pca/real_data/1KG_pruned_forPCA.gds")
geno <- index.gdsn(gfile, "genotype")

d <- 2504
N <- 168047
n <- round(N/100)
ni <- if (i < 99) n else (N - 99*n)
xi <- read.gdsn(geno)
ind <- apply(xi,2, function(x){any(x>2 & x <0)})
xi<-xi[,!ind]
xi <- scale(xi)
ts <- Sys.time()
xi <- xi[,(((i-1)*n + 1):( (i-1)*n + ni) )]
####Step 1: compute parallel sketches corresponding to each data split
set.seed(l)
omega <- matrix(rnorm(d*p),d , p)
Y <- t(xi)%*%omega
Y <- xi %*% Y/N
te <- Sys.time()
#shat <- 0.02

#ll <- read.gdsn(geno, start=c(1, 1), count=c(30, -1))
#sigc <- ll %*% t(ll)/ N
runtime <- as.numeric(te - ts, unit = "secs")
fname<-paste(c("/ncf/xlin_covid/Users/sshen/fast_pca/real_data/results_",args,".RData"),collapse = '_')
save(Y, runtime,file=fname)