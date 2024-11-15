args=(commandArgs(TRUE))

d <- as.numeric(args[1])
mc <- as.numeric(args[2])
rt <- as.numeric(args[3])


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
#d = 500,1000
#d = 2000

K = 3
m = 20
ni = 1000
n = ni*m


Delta2 = n^{2/3}


p = 4*K
p0 = 4*K
#rt = 0.2
L = round(rt*d/p)
#kp = K+1

it = 1

q = 7
lmode <- rt < 1

set.seed(152)
Theta = matrix(rnorm(n*K), n, K)
Theta <- Theta * (Delta2/2/n)^{0.5}


F_st <- matrix(0,d,K)
nk <- round(d/K)


for (k in 1:(K-1)){
  F_st[(((k-1)*nk+1):(k*nk)),k] <- 1
}
F_st[(((K-1)*nk+1):d),K] <- 1

svd1 <- svd(F_st)

mid <- t(Theta) %*% Theta 
svd2 <- svd(mid)

Vk <- svd1$u %*% svd2$u

EX <- Theta %*% t(F_st)


set.seed(mc+10086)

X <- matrix(rnorm(n*d),n,d) + EX

runtime <- 0




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
    y <- t(xi)%*%(xi%*%omega)
    Y <- Y+y
    te <- Sys.time()
    tm_ls <- c(tm_ls,difftime(te,ts,units="secs"))
  }
  tm <- max(tm_ls)
  ts <- Sys.time()
  Y <- Y-n*omega
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
    Bl <- t(yl_ls[[l]])%*%vtild/sqrt(p)/n
    svdl <- svd(Bl)
    Bl <- svdl$u %*% diag((svdl$d)^{-1}) %*% t(svdl$v)/n
    Bo <- rbind(Bo, Bl)
    Ysig <- cbind(Ysig, (n*omega_ls[[l]]+yl_ls[[l]])/sqrt(p))
    Omega <- cbind(Omega, omega_ls[[l]]/sqrt(p))
  }
  Sig_hat <- t(Omega)%*%Ysig
  Sig_hat <- t(Bo)%*%Sig_hat %*% Bo *(d*n^{1/3}/L)
  sc <- (L*d*n^{1/3})^{0.5}
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
  lmtilde <- lmtilde - n*diag(rep(1,K))
  lminv <- solve(lmtilde)
  Sig_hat <- (lminv+n*lminv %*% lminv)*(d^2*n^{1/3})
  sc <- d*n^{1/6}
  te <- Sys.time()
  tsig <- tsig + difftime(te,ts,units="secs")
}
runtime <- runtime + tsig

fname<-paste(c("/ncf/xlin_covid/Users/sshen/dissertation_1.5/example_3/results1/results_",args,".RData"),collapse = '_')
save(vtild,Sig_hat, sc, runtime,file=fname)


