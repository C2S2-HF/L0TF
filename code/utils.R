library(limSolve)
library(ggplot2)
library(ggpubr)

add_ind <- function(c1, c2){
  for(i in c1){
    if(i %in% c2) break
  }
  i
}

DiffMat1 <- function(n){
  D = cbind(-diag(n-1),0) + cbind(0, diag(n-1))
  return(D)
}
DiffMat <- function(n, k=1){
  D = DiffMat1(n)
  if (k>1){for (l in 1:(k-1)) D=DiffMat1(n-l)%*%D}
  return(D)
}


SimuEx <- function(n, sigma=0.1, q=0, nknot=4, seed=NA, RandKnot=FALSE, AdaKnot=FALSE) {
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  A=round(seq(0, n, length.out=nknot+2))[seq(2,nknot+1)]
  if(RandKnot) A = sort(sample(seq(6, n-5, 5), nknot))
  if(AdaKnot) A = round(seq(1, sqrt(n), length.out=nknot+2)^2)[seq(2,nknot+1)]
  tau = x[A]
  tau1 = c(0, tau, 1)
  if (q==0){
    aa = 1-seq(1,nknot+1)%%2
    y0 = 0*x
    for (j in 1:(nknot+1)) y0[x>tau1[j] & x<=tau1[j+1]] = aa[j]
  }
  if (q==1) {
    aa = 2*(-1)^seq(1,nknot+1)
    phi = rep(1, n)
    for (j in 1:(nknot+1)) phi = cbind(phi,pmin(pmax(x-tau1[j], 0), tau1[j+1]-tau1[j]))
    y0 = phi%*%c(0.5+1/(nknot+1),aa)
  }
  y = y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau, SetA = A, n=n, q=q, nknot=nknot, sigma=sigma))  
}


plt <- function(data, resL0){
  ltype = ifelse(data$q==0, "s", "l")
  plot(data$x, data$y, type='p', col='grey', xlab="", ylab="")
  lines(data$x, data$y0, col=1, lty=2, lwd=1, type=ltype) 
  lines(data$x, resL0$alpha, col=2, lwd=1, type=ltype)
  knot <- data$x[resL0$A+data$q]
  points(knot, rep(par("usr")[4], length(knot)), pch=3, cex=2, col=2)
  
}

inout_loop <- function(l1, title, xlab, ylab){
  df.m <- apply(l1, 2, mean)
  df.s <- apply(l1, 2, sd)
  dat <- data.frame(iter=1:length(df.m), df.m=df.m, df.s=df.s)

  plt <- ggplot(dat, aes(x=iter, y = df.m)) +
    geom_point(position = position_dodge(0)) +
    geom_line(aes(y = df.m)) +
    geom_errorbar(aes(ymin = df.m - df.s, ymax = df.m + df.s), 
                  width = 0.8, position = position_dodge(0))+
    ggtitle(title)+
    labs(x=xlab, y=ylab) +
    scale_x_continuous(breaks = seq(0, length(df.m), by = 1)) +
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.title.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y = element_text(size = 10,color="black"),
          plot.title = element_text(size = 10, hjust = 0.5))
  plt
}

convsteps <- function(l1, nlist, title, xlab, ylab){
  df.m <- apply(l1, 1, mean)
  df.s <- apply(l1, 1, sd)
  dat <- data.frame(iter=nlist, df.m=df.m, df.s=df.s)
  plt <- ggplot(dat, aes(x=iter, y = df.m)) +
    geom_point(position = position_dodge(0)) +
    geom_errorbar(aes(ymin = df.m - df.s, ymax = df.m + df.s), 
                  position = position_dodge(0))+
    ggtitle(title)+
    labs(x=xlab, y=ylab) +
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.title.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y = element_text(size = 10,color="black"),
          plot.title = element_text(size = 10, hjust = 0.5))
  plt
}

timesteps <- function(l1, nlist, title, xlab, ylab){
  list1 <- c()
  for(i in 1:length(nlist)){
    list1 <- c(list1, rep(nlist[i], 100))
  }
  list2 <- c()
  for(i in 1:dim(l1)[2]){
    list2 <- c(list2, l1[,i])
  }
  dat <- data.frame(n=list1, time=list2)
  plt <- ggplot(dat, aes(x=n, y = time)) +
    geom_point() +
    geom_smooth(color="black") +
    ggtitle(title)+
    labs(x=xlab, y=ylab) +
    theme(axis.text.x = element_text(size = 10,color="black"),
          axis.title.x = element_text(size = 10,color="black"),
          axis.text.y = element_text(size = 10,color="black"),
          axis.title.y = element_text(size = 10,color="black"),
          plot.title = element_text(size = 10, hjust = 0.5))
    plt
}
ToyEx <- function(n, sigma=0.1, q=0, seed=NA) {
  tau = 0.5
  if (!is.na(seed)) set.seed(seed)
  x = seq(1/n, 1,length.out = n)
  if (q==0){
    y0 = 0*x; y0[x>tau] = 1
  }
  if (q==1) {
    y0 = 2*(tau-x); y0[x>tau] = 2*(x[x>tau]-tau)
  }
  y = y0 + sigma*rnorm(n)
  return(list(y = y, x = x, y0 = y0, tau = tau))  
}
