nn=seq(256, 2048, 256)
load("SimuResultBlocks.RData")
png(paste0("SimuResultBlocks.png"),   pointsize = 25, width=1000, height=1600)
op <- par(mfrow=c(1,1), mar=c(4,2,3,2))
m <- matrix(c(1:6,7,7),nrow = 4, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.3, 0.3, 0.3, 0.2))
tmp=c("MSE", "nknot", "dP", "dS", "dH", "time")
tit <- c("MSE", "df", "dP", "dS", "dH", "log(Time)")
for (i in 1:6) {
  if(i==2){
    matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                      ResultPt[,tmp[i]], ResultWs[,tmp[i]], ResultNt[,tmp[i]], ResultLc[,tmp[i]]),
            type="b", lty=1, pch=c(1,2,3,4,5,6,7), col=c(2,4,3,5,6,7,8), ylim=c(5,39),
            main=tit[i], xlab="n", ylab="")
    abline(h=11, lty=2)
  }else{
    if(i!=6){
      matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                        ResultPt[,tmp[i]], ResultWs[,tmp[i]], ResultNt[,tmp[i]], ResultLc[,tmp[i]]),
              type="b", lty=1, pch=c(1,2,3,4,5,6,7), col=c(2,4,3,5,6,7,8),
              main=tit[i], xlab="n", ylab="")
    }else{
        matplot(nn, log10(cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                                ResultPt[,tmp[i]], ResultWs[,tmp[i]], ResultNt[,tmp[i]], ResultLc[,tmp[i]])),
                type="b", lty=1, pch=c(1,2,3,4,5,6,7), col=c(2,4,3,5,6,7,8),
                main=tit[i], xlab="n", ylab="")

      }
  }

  #legend("topright", c("L0TF", "L1TF", "Spline"), lty=1, pch=c(1,2,3), col=c(2,4,3))
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top", legend=c(expression(l[0]-TF), expression(l[1]-TF), "pspline", "PELT", "WBS", "NOT", "L0tfc"), lty=1, pch=c(1,2,3,4,5,6,7),
       cex=1.2, col=c(2,4,3,5,6,7,8), horiz = TRUE)
par(op)
# par(op)
dev.off()


load("SimuResultWave.RData")
png(paste0("SimuResultWave.png"),   pointsize = 25, width=1000, height=1600)
op <- par(mfrow=c(1,1), mar=c(4,2,3,2))
m <- matrix(c(1:6,7,7),nrow = 4, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.3, 0.3, 0.3, 0.2))
tmp=c("MSE", "nknot", "dP", "dS", "dH", "time")
tit <- c("MSE", "df", "dP", "dS", "dH", "log(Time)")
for (i in 1:6) {
  if(i==2){
    matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                      ResultCp[,tmp[i]], ResultNt[,tmp[i]]),
            type="b", lty=1, pch=c(1,2,3,4,5), col=c(2,4,3,1,5),
            main=tit[i], xlab="n", ylab="")
    abline(h=8, lty=2)
  }else{
    if(i!=6){
      matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                        ResultCp[,tmp[i]], ResultNt[,tmp[i]]),
              type="b", lty=1, pch=c(1,2,3,4,5), col=c(2,4,3,1,5),
              main=tit[i], xlab="n", ylab="")
    }else{
      matplot(nn, log10(cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                              ResultCp[,tmp[i]], ResultNt[,tmp[i]])),
              type="b", lty=1, pch=c(1,2,3,4,5), col=c(2,4,3,1,5),
              main=tit[i], xlab="n", ylab="")

    }
  }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top", legend=c(expression(l[0]-TF), expression(l[1]-TF), "pspline","CPOP","NOT"), lty=1, pch=c(1,2,3,4,5),
       cex=1.2, col=c(2,4,3,1,5), horiz = TRUE)
par(op)
# par(op)
dev.off()


load("SimuResultDoppler.RData")
png(paste0("SimuResultDoppler.png"),  pointsize = 22, width=1000, height=700)
op <- par(mfrow=c(1,1), mar=c(4,2,3,2))
m <- matrix(c(1:2,3,3),nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.7,0.3))
tmp <- c("MSE", "time")
tit <- c("MSE", "log(Time)")
i=1
matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]], ResultNt[,tmp[i]]),
        type="b", lty=1, pch=c(1,2,3,4), col=c(2,4,3,1),
        main=tit[i], xlab="n", ylab="")
i=2
matplot(nn, log10(cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]], ResultNt[,tmp[i]])),
        type="b", lty=1, pch=c(1,2,3,4), col=c(2,4,3,1),
        main=tit[i], xlab="n", ylab="")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "center", legend=c(expression(l[0]-TF), expression(l[1]-TF), "pspline", "NOT"), lty=1, pch=c(1,2,3,4),
       cex=1.2, col=c(2,4,3,1), horiz = TRUE)
par(op)
dev.off()
