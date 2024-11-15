args=(commandArgs(TRUE))
mc <- as.numeric(args[1])
K = 4
p = 50
m = 10
load("/ncf/xlin_covid/Users/sshen/dissertation_1.5/real_data/1000g_sbm95.RData")
info <- read.delim("/n/home_fasse/shutingshen/ncf_user/fast_pca/1KG_TRACE_pca.txt", sep = " ")

ll <- order(info$Population.2)
pop_vec <- info$Population.2[ll]

loc_pure <- pop_vec %in% c("AFR", "EAS", "EUR", "SAS")

Mh <- sbm95[loc_pure, loc_pure]
d = dim(Mh)[1]
dj <- d/m
thetah = sum(Mh)/d^2/2
mu = d*log(d)*sqrt(thetah/p)/12

t_ls <- c()
Yt <- 0
for (j in 1:m){
  loc <- (dj*(j-1)+1):(dj*j)
  ts <- Sys.time()
  omega <- matrix(rnorm(dj*p),dj,p)
  Ytj <- Mh[,loc]%*%omega
  Yt <- Yt + Ytj
  te <- Sys.time()
  t_ls <- c(t_ls, difftime(te,ts,units = "secs"))
}

ts <- Sys.time()
svdY <- svd(Yt/sqrt(p))
te <- Sys.time()
rt <- difftime(te,ts,units = "secs") + max(t_ls)
vkl <- svdY$u[,1:K] 
#####estimate K
dd <- svdY$d
diff <- dd-dd[p]
k <- 1
flag <- FALSE
while (!flag) {
  flag <- all(diff[k:p] <= mu)
  k<-k+1
}
kl<-k-2



fname<-paste(c("/ncf/xlin_covid/Users/sshen/dissertation_1.5/real_data/results/results_",args,".RData"),collapse = '_')

save(vkl,kl,rt, file = fname)