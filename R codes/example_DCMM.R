args=(commandArgs(TRUE))
d <- as.numeric(args[1])
mc <- as.numeric(args[2])
rt <- as.numeric(args[3])

fast_pca_final<-function(y_ls,p0,qq){
  d < -dim(y_ls[[1]])[1]
  omega<-matrix(rnorm(d*p0),d,p0)
  
  for (q in 1:qq){
    Yt <- 0
    L <- length(y_ls)
    for (l in 1:L){
      Y <- t(y_ls[[l]])%*%omega
      Y <- y_ls[[l]] %*% Y
      Yt <- Yt + Y
    }
    Yt <- Yt / L
    omega <- Yt
  }
  
  
  Q<-svd(Yt)$u
  return(list('u'=Q))
}


K =3
a = 0.7
b = 0.2


p = 4*K
p0 = 4*K
#rt = 10
L = round(rt*d/p)
lmode = 1
q = 7
it = 1



ext = 12

set.seed(153)
Pi <- matrix(runif(K*K),K,K)

Pi <- Pi + diag(rep(10,K))

Pi <- Pi / rowSums(Pi)

Pi1 <- matrix(0,d,K)
nk <- round(d/K)


for (k in 1:(K-1)){
  Pi1[(((k-1)*nk+1):((k-1)*nk+round(nk/3))),] <- rep(Pi[k,],each = round(nk/3))
  Pi1[(((k-1)*nk+round(nk/3)+ 1):(k*nk)),k] <- 1
  
}
Pi1[(((K-1)*nk+1):((K-1)*nk+round(nk/3))),] <- rep(Pi[K,],each = round(nk/3))
Pi1[(((K-1)*nk+round(nk/3)+1):d),K] <- 1


P <- diag(rep(0.3,3)) + 0.1*c(2,1,1)%*%t(c(2,1,1)) + 0.1 * c(1,-1,-1)%*%t(c(1,-1,-1))
svd1 <- svd(Pi1)

mid <- diag(svd1$d) %*% t(svd1$v) %*% P %*% svd1$v %*% diag(svd1$d)
svd2 <- svd(mid)

Vk <- svd1$u %*% svd2$u


M <- Pi1 %*% P %*% t(Pi1)
set.seed(mc)
Mh <- matrix(rbinom(rep(1,d^2),1,c(M)),d,d)
Mh[lower.tri(Mh)] <- t(Mh)[lower.tri(Mh)]




runtime <- 0

y_ls <- list()
t_ls <- c()

for (l in 1:L){
  ts <- Sys.time()
  omega <- matrix(rnorm(d*p),d,p)
  Y <- Mh%*%omega
  vl <- svd(Y)$u[,1:K]
  y_ls <- append(y_ls,list(vl))
  te <- Sys.time()
  t_ls <- c(t_ls,  difftime(te,ts,units="secs"))
}

runtime <- runtime + max(t_ls)

ts <- Sys.time()
vtild <- fast_pca_final(y_ls, p0,q)$u[,1:K]

if (lmode == 1){
  
  
  lmtild <- (t(vtild)%*%Mh%*%vtild)
  Mtild <- vtild%*%lmtild%*%t(vtild)
  
  Sig_hat1 <- solve(lmtild)%*%t(vtild)%*%( (    (Mtild[1,])*(1-Mtild[1,])   ) *vtild  )%*%solve(lmtild)
  Sig_hat2 <- solve(lmtild)%*%t(vtild)%*%( (   (Mtild[2])*(1-Mtild[2])   ) *vtild  )%*%solve(lmtild)
  Sig_hat <- Sig_hat1 + Sig_hat2
  }
te <- Sys.time()

runtime <- runtime +  difftime(te,ts,units="secs")

vkh <- svd(Mh)$u[,1:K]

fname<-paste(c("/ncf/xlin_covid/Users/sshen/dissertation_1.5/example_2/results/results_",args,".RData"),collapse = '_')
save(vtild,Sig_hat,Sig_hat1,Sig_hat2, vkh,runtime,file=fname)
