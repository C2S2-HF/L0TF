##################################################################################################
###                                                                                            ###
### Replication code for "L0 trend filtering" by Canhong Wen, Xueqin Wang and Aijun Zhang      ###
### This file generates Figures B.10, B.12 and B.14 in Appendix B.2.                           ###
### Updated on 10 April 2023                                                                   ###
###                                                                                            ###
##################################################################################################


## Figure B.10---------------------------
nn=seq(32,256,32)
mode = "pre"
load(paste0("plot_results/",mode,"_Blocks.RData"))
png(paste0("plots/",mode,"_Blocks.png"),   pointsize = 25, width=900, height=1100)
op <- par(mfrow=c(1,1), mar=c(4,2,3,2))
m <- matrix(c(1:6,7,7),nrow = 4, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.3, 0.3, 0.3, 0.2))
tmp=c("MSE", "nknot", "dP", "dS", "dH", "time")
tit <- c("MSE", "df", "dP", "dS", "dH", "log(Time)")
for (i in 1:6) {
  if(i==2){
    matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                      ResultPt[,tmp[i]], ResultWs[,tmp[i]], ResultNt[,tmp[i]], ResultLc[,tmp[i]]),
            type="b", lty=1, pch=c(1,2,3,4,5,6,7), col=c(2,4,3,5,6,7,8), ylim=c(1,39),
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
                                ResultPt[,tmp[i]], ResultWs[,tmp[i]], 
                                ResultNt[,tmp[i]], ResultLc[,tmp[i]])),
                type="b", lty=1, pch=c(1,2,3,4,5,6,7), col=c(2,4,3,5,6,7,8),
                main=tit[i], xlab="n", ylab="")

      }
  }
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top", legend=c(expression(l[0]-TF), expression(l[1]-TF),
                           "pspline", "PELT", "WBS", "NOT", expression(l[0]-MIP)),
       lty=1, pch=c(1,2,3,4,5,6,7), cex=1, col=c(2,4,3,5,6,7,8), horiz = TRUE)
par(op)
dev.off()

## Figure B.12---------------------------
nn=seq(32,256,32)
load(paste0("plot_results/",mode,"_Wave.RData"))
rm(ResultCp)
png(paste0("plots/",mode,"_Wave.png"),   pointsize = 25, width=900, height=1100)
op <- par(mfrow=c(1,1), mar=c(4,2,3,2))
m <- matrix(c(1:6,7,7),nrow = 4, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.3, 0.3, 0.3, 0.2))
tmp=c("MSE", "nknot", "dP", "dS", "dH", "time")
tit <- c("MSE", "df", "dP", "dS", "dH", "log(Time)")
for (i in 1:6) {
  if(i==2){
    matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                      ResultNt[,tmp[i]], ResultLc[,tmp[i]]),
            type="b", lty=1, pch=c(1,2,3,6,7), col=c(2,4,3,7,8),
            main=tit[i], xlab="n", ylab="")
    abline(h=8, lty=2)
  }else{
    if(i!=6){
      matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                        ResultNt[,tmp[i]], ResultLc[,tmp[i]]),
              type="b", lty=1, pch=c(1,2,3,6,7), col=c(2,4,3,7,8),
              main=tit[i], xlab="n", ylab="")
    }else{
      matplot(nn, log10(cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]], ResultSp[,tmp[i]],
                              ResultNt[,tmp[i]], ResultLc[,tmp[i]])),
              type="b", lty=1, pch=c(1,2,3,6,7), col=c(2,4,3,7,8),
              main=tit[i], xlab="n", ylab="")

    }
  }
}
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top", legend=c(expression(l[0]-TF), expression(l[1]-TF), "pspline","NOT",expression(l[0]-MIP)), 
       lty=1, pch=c(1,2,3,6,7), cex=1, col=c(2,4,3,7,8), horiz = TRUE)
par(op)
dev.off()

## Figure B.14---------------------------
load(paste0("plot_results/",mode,"_Doppler.RData"))
png(paste0("plots/",mode,"_Doppler.png"),  pointsize = 22, width=1000, height=600)
op <- par(mfrow=c(1,1), mar=c(4,2,3,2))
m <- matrix(c(1:2,3,3),nrow = 2, ncol = 2, byrow = TRUE)
layout(mat = m, heights = c(0.6,0.3))
tmp <- c("MSE", "time")
tit <- c("MSE", "log(Time)")
i=1
matplot(nn, cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]],
                  ResultSp[,tmp[i]], ResultNt[,tmp[i]], ResultLc[,tmp[i]]),
        type="b", lty=1, pch=c(1,2,3,6,7), col=c(2,4,3,7,8),
        main=tit[i], xlab="n", ylab="")
i=2
matplot(nn, log10(cbind(ResultL0[,tmp[i]], ResultL1[,tmp[i]],
                        ResultSp[,tmp[i]], ResultNt[,tmp[i]], ResultLc[,tmp[i]])),
        type="b", lty=1, pch=c(1,2,3,6,7), col=c(2,4,3,7,8),
        main=tit[i], xlab="n", ylab="")
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top", legend=c(expression(l[0]-TF), expression(l[1]-TF), "pspline", "NOT", expression(l[0]-MIP)), 
       lty=1, pch=c(1,2,3,6,7), cex=1, col=c(2,4,3,7,8), horiz = TRUE)
par(op)
dev.off()

