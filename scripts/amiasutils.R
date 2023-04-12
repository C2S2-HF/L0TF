library(limSolve)

amias_R <- function(y, D, k, A, rho, q){
  n <- length(y)
  set <- nrow(D)
  S <- 1:set
  I <- setdiff(S, A)
  for(m in 1:20){
    preA <- A
    u <- rep(0, set)
    eps <- rep(0, set)

    if(!is.null(I)){
      DI2 <- D[I,]%*%t(D[I,])
      nup <- sum(DI2[1,-1]!=0)
      u[I] <- Solve.banded(DI2, nup, nup, B = D[I,]%*%y)
    }
    if(!is.null(A)) eps[A] <- D[A,]%*%(y-t(D)%*%u)

    A <- sort(order(abs(eps+u/rho),decreasing=TRUE)[1:k])
    I <- sort(setdiff(S, A))
    if(setequal(preA, A)) break
  }
  # if(m==20) warning("reach max iteration")
  
  alpha <- y - t(D[I,])%*%u[I]
  list(A=A, u=u, v=eps, alpha=alpha, iter=m)
}

samias_R <- function(y, D, kmax, rho, q, eps=0.001){
  n <- length(y)
  u <- solve(D%*%t(D))%*%D%*%y
  set <- n-q-1
  A <- c()
  S <- 1:set
  I <- setdiff(S, A)
  bic <- c()
  
  docA <- list()
  docv <- list()
  docu <- list()
  docalpha <- list()
  dociter <- c()
  docmse <- c()
  for(k in 1:kmax){
    ord_u <- order(abs(u), decreasing=TRUE)
    A <- sort(union(A, add_ind(ord_u, I)))
    inloop <- amias_R(y, D, k, A, rho, q)
    A <- inloop$A
    u <- inloop$u
    mse <- sum((y-inloop$alpha)^2)/n
    bic <- c(bic, n*log(mse)+2*(q+1+k)*log(n))
    docA[[k]] <- A
    docv[[k]] <- inloop$v
    docu[[k]] <- inloop$u
    docalpha[[k]] <- inloop$alpha
    dociter <- c(dociter, inloop$iter)
    docmse <- c(docmse, mse)
    if(mse<=eps) break
  }
  idx <- which.min(bic)
  list(A=docA[[idx]], alpha=docalpha[[idx]], v=docv[[idx]], u=docu[[idx]], knot=idx, iters=dociter, mse=docmse, outiter=k)
}


ramias <- function(data, kmax){
    iters <- c()
    mse <- c()
    for(k in 1:kmax){
        setA <- sort(sample(seq(6, data$n-5, 5), k))
        res = amias(y = data$y, D = DiffMat(data$n, data$nknot), A = setA, k = k, rho = data$n**(data$q+1), q = data$q)
        iters <- c(iters, res$iter)
        mse <- c(mse, sum((data$y-res$alpha)^2)/data$n)
    }
    list(mse=mse, iters=iters)
}
