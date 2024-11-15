args=(commandArgs(TRUE))
##args=c(omega_l)
## m = 100, L = 80, p = 50 

l <- as.numeric(args[1])
p <- as.numeric(args[2])

d <- 2504
N <- 168047
mu= (d/(N*p)^{.5}*log(d))^{0.75}/12
###### Step 2 aggregating sketches across distributed splits
Yt <- 0 
rt <- c()
for (i in 1:100) {
  fname <- paste(c("/ncf/xlin_covid/Users/sshen/fast_pca/real_data/results_",c(i,l,50),".RData"),collapse = '_')
  load(fname)
  Yt <- Yt+Y
  rt <- c(rt, runtime)
}

rt <- max(rt)
ts<- Sys.time()
set.seed(l)
omega <- matrix(rnorm(d*p),d , p)

Yt <- Yt - 0.730152*omega

##### Compute parallel PCA
vk_hat <- prcomp(Yt)$x[,1:25]
te <- Sys.time()

rt <- rt + te - ts

##### Local estimation of K
dd <- svd(Yt/sqrt(p))$d
diff <- dd-dd[p]
k <- 1
flag <- FALSE
while (!flag) {
  flag <- all(diff[k:p] <= mu)
  k<-k+1
}
kl<-k-2

fname<-paste(c("/ncf/xlin_covid/Users/sshen/fast_pca/real_data/results_collect1",args,".RData"),collapse = '_')
save(vk_hat,rt,kl, runtime,file=fname)
