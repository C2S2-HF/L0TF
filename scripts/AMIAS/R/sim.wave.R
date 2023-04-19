# simulate wave
# n determines length
# sigma is residual sd
# output is data (y) and mean function (x) and changepoints (cpt)
sim.wave <- function(n, sigma = 1){
  
  x.breaks <- c(0, 256, 512, 768, 1024, 1152, 1280, 1344, 1408, 1472)/1472*n
  # x.breaks <- c(0, 256, 512, 768, 1024, 1152, 1280, 1344,  1500)/1500*n
  m <- length(x.breaks)-2
  x.breaks <- round(x.breaks)
  x.breaks[1] <- 1
  h <- (-1)^(1:m) *(1:m)/512
  x <- rep(0,n)
  slope <- 1/512 # wave1: 1/256 # wave2: 1/512 wave3: 1/384
  x[x.breaks[1]:x.breaks[2]] <- 1 + x.breaks[1]:x.breaks[2] * slope
  
  for(j in 2:(m+1)) {
    slope <- slope +  h[j-1]
    for(kk in x.breaks[j]:x.breaks[j+1]) x[kk] <- x[kk-1] + slope
  }
  
  a <- 3*sigma/sd(x)  # the signal to noise ratio is 3.
  x <- a*x
  y <- x + sigma*rnorm(n)
  
  cpt <- x.breaks[2:(m+1)]
  
  
  
  return(list(y = y, y0 = x, knots = cpt))
}
