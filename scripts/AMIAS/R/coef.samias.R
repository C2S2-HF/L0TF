coef.samias <- function(object,type = c("coef","primal", "dual", "active"), k, ...) {
  type <- match.arg(type)
  if (!(type %in% c("coef","primal", "dual", "active"))) {
    stop("Invalid type, must be \"coef\", or \"primal\", or \"dual\", or \"active\".")
  }
  
  if(missing(k)) {
    # Nothing was specified, use optimal k
    if (type == "coef") return(list(alpha=as.matrix(object$alpha), k=object$kopt)) 
    if (type == "primal") return(list(v=as.matrix(object$v),k=object$kopt))
    if (type == "dual") return(list(u=as.matrix(object$u),k=object$kopt))
    if (type == "active") return(list(A=object$A, k=object$kopt))
  }
  else{
    if (type == "coef") return(list(beta=object$alpha.all[,k], k=k)) 
    if (type == "primal") return(list(v=object$v.all[,k],k=k))
    if (type == "dual") return(list(u=object$u.all[,k],k=k))
    if (type == "active") return(list(A=object$A.all[k], k=k))
  }

}

