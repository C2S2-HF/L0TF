library(gurobi)
library(pracma)
SimuBlocks <- function (n, sigma = 0.1, seed=NA) 
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 
           0.76, 0.78, 0.81)
  h <- c(0, 4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)/5
  A = sapply(tau, function(z) which(x>=z)[1])
  tau1=c(0,tau,1)
  y0 = 0*x
  for (j in 1:(length(A)+1)) y0[x>tau1[j] & x<=tau1[j+1]] = h[j]
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))  
}

SimuWave <- function (n, sigma = 0.1, seed=NA) 
{
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  tau = c(256, 512, 768, 1024, 1152, 1280, 1344, 1408)/1472
  nknot = length(tau)
  A = sapply(tau, function(z) which(x>=z)[1])
  tau=x[A]
  tau1 = c(0, tau, 1)
  h =  cumsum(c(1, (-1)^(1:nknot)*(1:nknot)))
  phi = rep(1, n)
  for (j in 1:(nknot+1)) phi = cbind(phi,pmin(pmax(x-tau1[j], 0), tau1[j+1]-tau1[j]))
  y0 = as.vector(phi%*%c(0, h))
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A))  
}

SimuDoppler <- function(n, sigma = 0.1, seed=NA){
  if (!is.na(seed)) set.seed(seed)
  x <- seq(1/n, 1,length.out = n)
  y0 <- sqrt(x*(1-x))*sin(2*pi*(1+0.05)/(x+0.05))
  y <- y0 + sigma*rnorm(n)
  return(list(y = y, y0 = y0, x=x))
}

# -------------------------------------
# Evalulate metrics for benchmark comparison
# beta is the estimate 
# y0 is the true signal
# cpts is the estimated change points
# tcpt is the true change points
# -------------------------------------
EvalMetrics <- function(beta, cpts=NULL, y0, tcpt = NULL){
  mse <- mean((beta-y0)^2)
  mad <- mean(abs(beta-y0))
  if(is.null(tcpt)){
    tab <- data.frame(MSE=mse, MAD=mad)
  }else{
    n <- length(beta)
    n.cpts <- length(cpts)
    segments.endpoints.true <- sort(unique(tcpt))
    segments.endpoints.est <- sort(unique(cpts))
    distm <- abs(matrix(rep(segments.endpoints.est, length(segments.endpoints.true)), nrow=length(segments.endpoints.est))
                 -matrix(rep(segments.endpoints.true, length(segments.endpoints.est)), nrow=length(segments.endpoints.est), byrow=TRUE))
    screening.dist <- max(apply(distm, 2, min)) * 100 / n # min distance for each true cpt
    precision.dist <- max(apply(distm, 1, min)) * 100 / n # min distance for each estimated cpt
    haus.dist <- max(screening.dist, precision.dist)
    tab <- data.frame(MSE=mse, MAD=mad, dS=screening.dist, dP=precision.dist, dH=haus.dist, nknot=n.cpts)
    #names(tab) <- c("MSE", "MAD", "d_S", "d_P", "d_H", "n.cpts")
  }
  return(tab)
}

DiffMat1 <- function(n){
  D = cbind(-diag(n-1),0) + cbind(0, diag(n-1))
  return(D)
}
DiffMat <- function(n, k=1){
  D = DiffMat1(n)
  if (k>1){for (l in 1:(k-1)) D=DiffMat1(n-l)%*%D}
  return(D)
}

l0tfc <- function(y=y, D_matrix, l0constraint=l0constraint, M=M, l2penalty=l2penalty,...){
  
  n <- length(y)
  rD <- nrow(D_matrix)
  A1 <- cbind(D_matrix,-M*diag(rD))
  A2 <- cbind(-D_matrix,-M*diag(rD))
  A3 <- matrix(c(rep(0,n),rep(1,(rD))),nrow=1)
  A_matrix <- rbind(A1,A2,A3)
  A_c <- ncol(A_matrix)
  A_r <- nrow(A_matrix)
  
  b_vector <- c(rep(0,(A_r-1)),l0constraint)
  
  c_vector <- c(-2*y,rep(0,rD))
  
  Q_matrix <- matrix(0,nrow=A_c,ncol=A_c)
  Q_matrix[1:n,1:n] <- diag(n) + l2penalty*(t(D_matrix)%*%D_matrix)
  
  x_ub <- c(rep(max(y),n),rep(1,rD))
  x_lb <- c(rep(min(y),n),rep(0,rD))
  
  x_type <- rep('C',n)
  z_type <- rep('B',rD)
  
  model <- list()
  
  model$A          <- A_matrix
  model$rhs        <-  b_vector
  model$sense      <- c(rep('<',A_r-1),"<=")
  model$obj        <- c_vector
  model$Q          <- Q_matrix
  model$modelsense <- 'min'
  model$lb         <- x_lb
  model$ub         <- x_ub
  model$vtype      <- c(x_type,z_type)
  
  params <- list(OutputFlag=0)
  
  result <- gurobi(model, params)
  
  return(result)
  
  rm(model, result, params)
}

BIC_l0tfc <- function(y, D, kmax, q){
  bwl=3
  n = length(y)
  Y = movavg(y, bwl, type="s")
  M = max(abs(diff(x = Y, differences = 2)))
  bic <- c()
  fine <- matrix(0,nrow=n,ncol=kmax)
  knots <- list()
  for(k in 1:kmax){
    temp <- l0tfc(y,D,l0constraint=k,M,l2penalty=0)
    fine[,k] <- temp$x[1:n]
    knots[[k]] <- which(temp$x[(n+1):(2*n-1)]==1)+1
    sighat <- median(abs(diff(y, diff = q+1)))/(qnorm(3/4)*sqrt(choose(2*(q+1), q+1)))
    bbb <- sum((y-fine[,k])^2)/sighat^2 + 2*log(n)*(length(knots[[k]])+q+1)
    bic <- c(bic, bbb)
  }
  idx <- which.min(bic)
  list(knot=knots[[idx]], fit=fine[,idx], bic=bic[idx])
}

SingleRunL0tfcTF <- function(dgm = "Blocks", n=300, sigma=0.1, seed=NA){
  
  if (dgm=="Blocks"){
    data = SimuBlocks(n=n, sigma=sigma, seed=seed)
    q=0; kmax=20;
    tic = Sys.time()
    D_matrix = DiffMat1(n)
    resL0tfc = BIC_l0tfc(data$y,D_matrix,kmax,q)
    toc = Sys.time()
  }
  if (dgm=="Wave"){
    data = SimuWave(n=n, sigma=sigma, seed=seed)
    q=1; kmax=20;
    tic = Sys.time()
    D_matrix = DiffMat(n,q)
    resL0tfc = BIC_l0tfc(data$y,D_matrix,kmax,q)
    toc = Sys.time()
  }
  if (dgm=="Doppler"){
    data = SimuDoppler(n=n, sigma=sigma, seed=seed)
    q=2; kmax=20;
    tic = Sys.time()
    D_matrix = DiffMat(n,q) 
    resL0tfc = BIC_l0tfc(data$y,D_matrix,kmax,q)
    toc = Sys.time()
  }
  metric = EvalMetrics(resL0tfc$fit, resL0tfc$knot/n, data$y0, data$tau)
  metric$time = as.numeric(toc-tic,  units = "secs")
  return(metric)
}

# -------------------------------------
# Blocks - Monte Carlo Replicates
# -------------------------------------
library(foreach)
library(doMC)
registerDoMC(cores=5)

# nn=seq(32,256,32); sigma=0.1;
# ResultLc <- NULL
# for (n in nn){
#   cat("\nRunning n =", n)
#   TabLc = foreach(mcseed = 1:10, .combine = rbind) %dopar%
#     SingleRunL0tfcTF(dgm = "Blocks", n=n, sigma=sigma, seed=mcseed)
#   ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
# }
# save(ResultLc, file="lcBlocks.RData")
# 
# 
# nn=seq(32,160,32); sigma=0.1;
# ResultLc <- NULL
# for (n in nn){
#   cat("\nRunning n =", n)
#   TabLc = foreach(mcseed = 1:10, .combine = rbind) %dopar%
#     SingleRunL0tfcTF(dgm = "Wave", n=n, sigma=sigma, seed=mcseed)
#   ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
# }
# save(ResultLc, file="lcWave.RData")

nn=seq(32,96,32); sigma=0.1;
ResultLc <- NULL
for (n in nn){
  cat("\nRunning n =", n)
  TabLc = foreach(mcseed = 1:10, .combine = rbind) %dopar%
    SingleRunL0tfcTF(dgm = "Doppler", n=n, sigma=sigma, seed=mcseed)
  ResultLc = rbind(ResultLc, c(n=n, colMeans(TabLc)))
}

save(ResultLc, file="lcDoppler.RData")
