print.samias <- function (object, digits = max(3, getOption("digits") - 3), ...) 
{
  cat("\nCall: ", deparse(object$call), "\n\n")
  print(cbind(k = object$k, Df = object$df, `BIC` = signif(object$bic, digits)))
}