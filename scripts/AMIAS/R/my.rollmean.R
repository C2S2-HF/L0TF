my.rollmean=function(y, h = 5, ...){
  z.m <- zoo::rollmean(y,k=h*2+1,...)
  nei.l <- cumsum(y[1:(2*h)])[(h+1):(2*h)]
  z.l <- nei.l / ((h+1):(2*h))
  y1 <- rev(y)
  nei.r <- cumsum(y1[1:(2*h)])[(h+1):(2*h)]
  z.r <- nei.r / ((h+1):(2*h))
  z.r <- rev(z.r)
  z <- c(z.l, z.m, z.r)
  z
}