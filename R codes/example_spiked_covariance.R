args=(commandArgs(TRUE))

d <- as.numeric(args[1])
mc <- as.numeric(args[2])
kp <- as.numeric(args[3])
rt <- as.numeric(args[4])
fast_pca_final<-function(y_ls,p0,qq){
  d <- dim(y_ls[[1]])[1]
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


K = 3
m = 20
ni = 1000
n = ni*m

sigma2 = 1
det = 2
p = 4*K
p0 = 4*K
L = round(rt*d/p)

it = 1

q = 7
lmode <- rt < 1

lambda = det*(K:1)

set.seed(151)
u <- matrix(rnorm(d*K), d,K)
u <- svd(u)$u

Sig <- u %*% diag(lambda) %*% t(u) + sigma2*diag(rep(1,d))

Sigsqrt <- svd(Sig)$u %*% diag(c((lambda+1)^{.5},rep(1,d-K)))


set.seed(mc+10086)

X <- matrix(rnorm(n*d),n,d)%*%t(Sigsqrt)

runtime <- 0

######estimate of sig2 and u_{K+1}
ts <- Sys.time()
SigK <- t(X[,1:kp]) %*% X[,1:kp]/n
svdk <- svd(SigK)
sh <- svdk$d[kp]
te <- Sys.time()
runtime <- runtime + difftime(te,ts,units="secs")

omega_ls <- list()
yl_ls <- list()
y_ls <- list()

tL_ls <- c()

###########distributed fast sketches
for (l in 1:L){
  omega <- matrix(rnorm(d*p),d,p)
  omega_ls <- append(omega_ls,list(omega))
  Y<- 0 
  tm_ls <- c()
  for(i in 1:m){
    ts <- Sys.time()
    xi <- X[(((i-1)*ni+1):(i*ni)),]
    y <- t(xi)%*%(xi%*%omega)/n
    Y <- Y+y
    te <- Sys.time()
    tm_ls <- c(tm_ls,difftime(te,ts,units="secs"))
  }
  tm <- max(tm_ls)
  ts <- Sys.time()
  Y <- Y-sh*omega
  yl_ls <- append(yl_ls,list(Y))
  vl <- svd(Y)$u[,1:K]
  y_ls <- append(y_ls,list(vl))
  te <- Sys.time()
  tL_ls <- c(tL_ls, tm + difftime(te,ts,units="secs"))
}

tL <- max(tL_ls)
runtime <- runtime + tL

ts <- Sys.time()
vtild <- fast_pca_final(y_ls, p0,q)$u[,1:K]
te <- Sys.time()
runtime <- runtime + difftime(te,ts,units="secs")



########estimation of covariance

if(lmode == 1){
  
  Bo <- c()
  Ysig <- c()
  Omega <- c()
  ts <- Sys.time()
  for (l in 1:L){
    Bl <- t(yl_ls[[l]])%*%vtild/sqrt(p)
    svdl <- svd(Bl)
    Bl <- svdl$u %*% diag((svdl$d)^{-1}) %*% t(svdl$v)
    Bo <- rbind(Bo, Bl)
    Ysig <- cbind(Ysig, sh*omega_ls[[l]]+yl_ls[[l]])
    Omega <- cbind(Omega, omega_ls[[l]])
  }
  Sig_hat <- t(Omega)%*%Ysig
  Sig_hat <- t(Bo)%*%Sig_hat %*% Bo /(d*L)
  sc <- (n*L*p/d)^{0.5}
  te <- Sys.time()
  tsig <- difftime(te,ts,units="secs")
}else{
  lmtilde <- 0
  tsig <- c()
  for (i in 1:m){
    ts <- Sys.time()
    xi <- X[(((i-1)*ni+1):(i*ni)),]
    lmtilde <- lmtilde + t(vtild)%*% t(xi) %*% (xi %*% vtild)
    te <- Sys.time()
    tsig <- c(tsig, difftime(te,ts,units="secs"))
  }  
  tsig <- max(tsig)
  
  ts <- Sys.time()
  lmtilde <- lmtilde/n - sh*diag(rep(1,K))
  lminv <- solve(lmtilde)
  Sig_hat <- sh*lminv+sh^2*lminv %*% lminv
  sc <- n^{.5}
  te <- Sys.time()
  tsig <- tsig + difftime(te,ts,units="secs")
}
runtime <- runtime + tsig

fname<-paste(c("/ncf/xlin_covid/Users/sshen/dissertation_1.5/example_1/results1/results_",args,".RData"),collapse = '_')
save(vtild,Sig_hat, sc, runtime,file=fname)

