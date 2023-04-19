coef.amias <- function(object, type=c("coef","primal", "dual", "active"), ...) {
  type <- match.arg(type)
  if (!(type %in% c("coef","primal", "dual", "active"))) {
    stop("Invalid type, must be \"coef\", or \"primal\", or \"dual\", or \"active\".")
  }
  
  if (type == "coef") return(list(alpha=as.matrix(object$alpha), k=object$k)) 
  if (type == "primal") return(list(v=as.matrix(object$v),k=object$k))
  if (type == "dual") return(list(u=as.matrix(object$u),k=object$k))
  if (type == "active") return(list(A=object$A, k=object$k))


}

