
# simulate blocks
# n determines length
# sigma is residual sd
# output is data (y) and mean function (x) and changepoints (cpt)
sim.blocks <- function(n, sigma = 1){
  x0 <- seq(0,1,length.out = n)
  x.breaks <- c(0.1, 0.13, 0.15,0.23,0.25,0.4,0.44, 0.65, 0.76, 0.78, 0.81)
  h <- c(4,-5,3,-4,5,-4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  x <- 0*x0
  for(i in 1:n){
    x[i] <- sum(h*(1 + sign(x0[i]-x.breaks))/2)    ## true signal
  }
  y <- x + sigma*rnorm(n)
  
  cpt <- 0*x.breaks
  for(i in 1:length(cpt)) 
    cpt[i] <- which.min(abs(x0-x.breaks[i]))
  
  return(list(y = y, y0 = x, knots = cpt))
}

