
samias <- function(y, D = NULL, D_type=c("tf0", "tfq", "user"), q = 0, kmax = min(20,m-1), rho = n^(q+1), tmax = 10, eps = 0, smooth = FALSE, h = 5, adjust = FALSE, delta = 10, adjust.max = 10, ...){
  
  D_type <- match.arg(D_type)
  # Checking
  if (storage.mode(y)=="integer") storage.mode(y) <- "double"
  if (any(is.na(y))) stop("Missing data (NA's) detected. Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing y to SAMIAS")
  n = as.integer(length(y))
  m = dim(D)[1]
  if(!is.null(D)){
    if(!(dim(D)[2]==n))stop("D and y do not have the same number of observations")
    if (storage.mode(D)=="integer") storage.mode(D) <- "double"
  }else{
    if(D_type == "user")stop("Choose a 'user' D_type without passing D to SAMIAS")
  }
  
  if(D_type != "tf0" && D_type != "tfq" && D_type != "user")stop("No defined method")
  
  
  # Whether to adjust the detected knots
  if(adjust) adjust = 1 else adjust = 0
  
  # Smooth
  if(smooth) y <- my.rollmean(y, h = h)
  
  
  kmax <- kmax + 1
  outer_itermax <- kmax
  tau <- 1
  if(is.null(D)){
    if(D_type == "tf0"){
      result = .Fortran("afusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),
                        matrix(0,nrow=n-1,ncol=outer_itermax),
                        matrix(0,nrow=n-1,ncol=outer_itermax),as.integer(tau),
                        as.integer(outer_itermax), as.integer(kmax),eps,rho,rep(as.integer(1),outer_itermax),as.integer(tmax),
                        as.integer(adjust), as.integer(adjust.max), as.integer(delta), PACKAGE="AMIAS")
      result[[3]] = result[[3]][,1:result[[7]]]
      result[[4]] = result[[4]][,1:result[[7]]]
      result[[5]] = result[[5]][,1:result[[7]]]
      result[[11]] = result[[11]][1:result[[7]]]
    }else if(D_type == "tfq"){
      if(q==0)stop("Parameter q hasn't set to gennerate D (polynomial trend filtering of order q)")
      q = as.integer(q)
      stopifnot(q>=1)
      vec = genDtf1d(k=q)
      veclen = as.integer(length(vec))
      DTD = .Fortran("dtdmul",vec,n,veclen,matrix(0,nrow=2*veclen-1,ncol=n-veclen+1))[[4]]
      result = .Fortran("batfusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),
                        matrix(0,nrow=n-veclen+1,ncol=outer_itermax),matrix(0,nrow=n-veclen+1,ncol=outer_itermax),
                        as.integer(tau),as.integer(outer_itermax),as.integer(kmax),eps,rho,rep(as.integer(1),outer_itermax),
                        as.integer(tmax),DTD,vec,veclen,
                        as.integer(adjust), as.integer(adjust.max), as.integer(delta), PACKAGE="AMIAS")
      result[[3]] = result[[3]][,1:result[[7]]]
      result[[4]] = result[[4]][,1:result[[7]]]
      result[[5]] = result[[5]][,1:result[[7]]]
      result[[11]] = result[[11]][1:result[[7]]]
    }
  }else{
    D_type = "user"
    DTD <- D%*%t(D)
    result = .Fortran("agfusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),matrix(0,nrow=m,ncol=outer_itermax),
                      matrix(0,nrow=m,ncol=outer_itermax),as.integer(tau),
                      as.integer(outer_itermax),as.integer(kmax),eps,rho,rep(as.integer(1),outer_itermax),
                      as.integer(tmax), DTD,D,as.integer(m),
                      as.integer(adjust), as.integer(adjust.max), as.integer(delta),PACKAGE="AMIAS")
    result[[3]] = result[[3]][,1:result[[7]]]
    result[[4]] = result[[4]][,1:result[[7]]]
    result[[5]] = result[[5]][,1:result[[7]]]
    result[[11]] = result[[11]][1:result[[7]]]
  }
  
  mse <- colMeans((y-result[[3]])^2)
  df = seq(from=tau,by=tau,length.out=result[[7]]) + q + 1
  bic <- n*log(mse) + 2*log(n)*df
  opt = which.min(bic)
  A.all = apply(result[[4]], 2, function(x) which(x!=0))
  
  output = structure(list(call = match.call(),
                          y = y,
                          D_type = D_type,
                          q = q,
                          k = seq(from=tau,by=tau,length.out=result[[7]]),
                          alpha = result[[3]][,opt,drop = TRUE],
                          v = result[[4]][,opt,drop = TRUE],
                          u = result[[5]][,opt,drop = TRUE],
                          A = A.all[[opt]],
                          df = df,
                          alpha.all = result[[3]],
                          v.all = result[[4]],
                          u.all = result[[5]],
                          A.all = A.all,
                          kopt = opt,
                          bic = bic,
                          iter = result[[11]],
                          smooth = smooth),
                     class = 'samias')
  
  output
  
}



