plot.amias <- function(object, add.knots = TRUE,...){
  
  n <- length(object$y)
  y0 <- seq(1/n, 1, length.out = n)
  yname = "y"
  if(object$smooth)yname <- "smooth-y"
  plot.default(y0, object$y, ylab = yname, col="lightgrey", pch = 20, xlab="",...)
  
  ## Plot the estimates
  type <- ifelse(object$q==0, 's', "l") # Specify the plot type for the estimate
  lines(y0, object$alpha, col = "red", type = type, lwd = 2, lty = 1)
  
  if(add.knots == TRUE){
    ## Add the locations for the detected knots
    knot <- y0[which(object$v!=0)]
    points(knot, rep(par("usr")[4], length(knot)), pch=3, cex=2, col="red")
  }
  
}