print.amias <- function(object, ...) {
  cat("\nCall:\n")
  dput(object$call)
  cat("\nOutput:\n")
  cat(paste("Generalized L0 coefficients with number of knots being ", object$k,".",
            "\n\n", sep=""))
}
