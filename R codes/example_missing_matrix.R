args=(commandArgs(TRUE))
d <- as.numeric(args[1])
mc <- as.numeric(args[2])
rt <- as.numeric(args[3])
m <- 10
dj <- d/m
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
det =2
lambda = det*(K:1)

p = 4*K
p0 = 4*K
L = round(rt*d/p)
lmode = (rt >= 1)
q = 7
it = 1

theta = 0.4
sig = d^{-1}*det*4



set.seed(154)

u <- matrix(rnorm(d*K), d,K)
u <- svd(u)$u

M <- u%*%diag(lambda)%*%t(u)


set.seed(mc)
E <- matrix(rnorm(d^2)*sig,d,d)
Mh <- M+E
obv_loc <- rbinom(d^2, 1, 1-theta) == 1
Mh[obv_loc] <- 0
Mh[lower.tri(Mh)] <- t(Mh)[lower.tri(Mh)]

runtime <- 0

yl_ls <- list()
y_ls <- list()
t_ls <- c()
omega_ls <- list()

for (l in 1:L){
tt_ls <- c()
  Y <- 0
omega <- c()
  for (j in 1:m){
    ts <- Sys.time()
    omegaj <- matrix(rnorm(d*p/m),d/m,p)   
 loc <- (dj*(j-1)+1):(dj*j)
    Yj <- Mh[,loc]%*%omegaj
    Y <- Y + Yj
    te <- Sys.time()
    tt_ls <- c(tt_ls,  difftime(te,ts,units="secs"))
omega <- rbind(omega,omegaj) 
 }

  ts <- Sys.time()
  vl <- svd(Y)$u[,1:K]
  te <- Sys.time()
  
  y_ls <- append(y_ls,list(vl))
  
  omega_ls <- append(omega_ls,list(omega))
  yl_ls <- append(yl_ls,list(Y))
  t_ls <- c(t_ls,  difftime(te,ts,units="secs") + max(tt_ls))

}

runtime <- runtime + max(t_ls)

ts <- Sys.time()
vtild <- fast_pca_final(y_ls, p0,q)$u[,1:K]
Ss <- sum(obv_loc & !lower.tri(Mh))
  that <- 1-2*Ss/d/(d+1)
  lmtild <- (t(vtild)%*%Mh%*%vtild)/that
  Mtild <- vtild%*%lmtild%*%t(vtild)
  sig2hat <- mean(((Mh - Mtild)[!obv_loc & !lower.tri(Mh)])^2)
  
if (lmode == 1){
  
  
  Sig_hat1 <- solve(lmtild)%*%t(vtild)%*%(((Mtild[1,]^2*(1-that) + sig2hat    )/that)*vtild)%*%solve(lmtild)
  Sig_hatd <- solve(lmtild)%*%t(vtild)%*%((( Mtild[d,]^2*(1-that) + sig2hat    )/that)*vtild)%*%solve(lmtild)
  Sig_hat2 <- solve(lmtild)%*%t(vtild)%*%((( Mtild[2,]^2*(1-that) + sig2hat    )/that)*vtild)%*%solve(lmtild)
  Sig_hat <- Sig_hat1 + Sig_hat2
  }else{
       Bo <- c()
  Omega <- c()

  for (l in 1:L){
    Bl <- t(yl_ls[[l]])%*%vtild/sqrt(p)
    svdl <- svd(Bl)
    Bl <- svdl$u %*% diag((svdl$d)^{-1}) %*% t(svdl$v)
    Bo <- rbind(Bo, Bl)
    Omega <- cbind(Omega, omega_ls[[l]]/sqrt(p))
  }
  
  Sig_hat1 <-  t(Bo)%*%t(Omega)%*%(((Mtild[1,]^2*(1-that) + sig2hat   )/that)*Omega)%*%Bo/L^2
  Sig_hatd <- t(Bo)%*%t(Omega)%*%((( Mtild[d,]^2*(1-that) + sig2hat    )/that)*Omega)%*%Bo/L^2
  Sig_hat2 <- t(Bo)%*%t(Omega)%*%((( Mtild[2,]^2*(1-that) + sig2hat    )/that)*Omega)%*%Bo/L^2
  Sig_hat <- Sig_hat1 + Sig_hat2
  
  }
te <- Sys.time()

runtime <- runtime +  difftime(te,ts,units="secs")

fname<-paste(c("/ncf/xlin_covid/Users/sshen/dissertation_1.5/example_4_new/results/results_",args,".RData"),collapse = '_')
if (lmode == 1){
save(vtild,lmtild,Mtild,Sig_hat, Sig_hat1, Sig_hatd, sig2hat,that, runtime,file=fname)
}else{
save(vtild,Bo,Omega,Mtild,Sig_hat, Sig_hat1, Sig_hatd, sig2hat,that, runtime,file=fname)
}
