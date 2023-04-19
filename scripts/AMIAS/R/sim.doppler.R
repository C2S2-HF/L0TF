# simulate doppler
# n determines length 
# sigma is residual sd
# output is data (y) and mean function (x) 
sim.doppler <- function(n, sigma = 1){
  x0 <- seq(0,1,length.out = n)
  x <- 10*sqrt(x0*(1-x0)) * sin(2*pi*(1+0.05)/(x0+0.05))
  y <- x + sigma*rnorm(n)
  return(list(y = y, y0 = x, knots = NULL))
}