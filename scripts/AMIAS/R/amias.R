amias <- function(y, D = NULL, D_type = c("tf0", "tfq", "user"), q = 0, k = 3, rho = n^(q+1), tmax = 10,
                  A = NULL, smooth = FALSE, h = 5, adjust = FALSE, delta = 10, adjust.max = 10,...){

  D_type <- match.arg(D_type)
  # Checking
  if (storage.mode(y)=="integer") storage.mode(y) <- "double"
  if (any(is.na(y))) stop("Missing data (NA's) detected. Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing y to AMIAS")
  n = as.integer(length(y))
  nm = dim(D)
  if(!is.null(D)){
    if(!(nm[2]==n))stop("D and y do not have the same number of observations")
    if (storage.mode(D)=="integer") storage.mode(D) <- "double"
  }else{
    if(D_type == "user")stop("Choose a 'user' D_type without passing D to AMIAS")
  }

  if(D_type != "tf0" && D_type != "tfq" && D_type != "user")stop("No defplot/ined method")

  # Whether to adjust the detected knots
  if(adjust) adjust = 1 else adjust = 0
  if(!is.null(A)&length(A)!=k) stop("The length of A should be either 0 or k!")
  if(is.null(A)){
    Anull = 1
    A = rep(0, k)
  } else {
    Anull = 0
  }
  # Whehter to smooth the data first
  if(smooth) y <- my.rollmean(y, h = h)

  DiffMat1 <- function(n){
    D = cbind(-diag(n-1),0) + cbind(0, diag(n-1))
    return(D)
  }
  DiffMat <- function(n, k=1){
    D = DiffMat1(n)
    if (k>1){for (l in 1:(k-1)) D=DiffMat1(n-l)%*%D}
    return(D)
  }
  # Whether the D matrix is given
  if(is.null(D)){
    if(D_type == "tf0"){
      rho = as.numeric(n)
      u = rep(0, n-1)
      I = setdiff(1:(n-1), A)

   #   DTD = .Fortran("dtdmul",vec,n,2,matrix(0,nrow=3,ncol=n-1),PACKAGE="AMIAS")[[4]]  # caculate D%*%t(D) and return in a banded form
      #if(Anull == 1){
        #D = DiffMat(n=n, k=1)
        #DTD <- D[I,]%*%t(D[I,])
        #nup <- sum(DTD[1,-1]!=0)
        #u[I] <- Solve.banded(DTD, nup, nup, B = D[I,]%*%y)
      #}
      result = .Fortran("fusedl0", n, y, rep(0,n), rep(0,n-1), u,
                        as.integer(A), as.integer(I), as.integer(Anull),
                        as.integer(k), rho, as.integer(1),as.integer(tmax),
                        as.integer(adjust), as.integer(adjust.max), as.integer(delta), PACKAGE="AMIAS")
    }else if(D_type == "tfq"){
      if(q==0)stop("Parameter q hasn't set to gennerate D (polynomial trend filtering of order q)")
      q = as.integer(q)
      stopifnot(q>=1)
      if(q==0) vec = c(-1, 1) else vec = drop(genDtf1d(k=q))
      veclen = as.integer(length(vec))
      u = rep(0, n-veclen+1)
      I = setdiff(1:(n-veclen+1), A)
      DTD = .Fortran("dtdmul",vec,n,veclen,matrix(0,nrow=2*veclen-1,ncol=n-veclen+1),PACKAGE="AMIAS")[[4]]  # caculate D%*%t(D) and return in a banded form
      #if(Anull == 1){
        #D = DiffMat(n=n, k=q+1)
        #DI2 <- D%*%t(D)
        #nup <- sum(DI2[1,-1]!=0)
        #u[I] <- Solve.banded(DI2, nup, nup, B = D[I,]%*%y)
     # }

      result = .Fortran("btfusedl0", n, y, rep(0,n), rep(0,n-veclen+1), u,
                        as.integer(A), as.integer(I), as.integer(Anull),
                        as.integer(k), rho, as.integer(1), as.integer(tmax), DTD, vec, veclen,
                        as.integer(adjust), as.integer(adjust.max), as.integer(delta), PACKAGE="AMIAS")
    }
  }else{
    D_type = "user"
    I = setdiff(1:nm[1], A)
    u = rep(0, nm[1])
    DTD = D%*%t(D)
    result = .Fortran("gfusedl0", n, y, rep(0,n), rep(0,nm[1]), u,
                      as.integer(A), as.integer(I), as.integer(Anull), as.integer(k),
                      rho, as.integer(1), as.integer(tmax), DTD, D, as.integer(nm[1]),
                      as.integer(adjust), as.integer(adjust.max), as.integer(delta), PACKAGE="AMIAS")
  }


  output = structure(list(call = match.call(),
                          y = y,
                          D_type = D_type,
                          q = q,
                          k = k,
                          alpha = result[[3]],
                          v = result[[4]],
                          u = result[[5]],
                          A = which(result[[4]]!=0),
                          df = k + q + 1,
                          iter = result[[11]],
                          smooth = smooth),
                     class = 'amias')

  output
}

