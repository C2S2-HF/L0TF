##################################################################################################
###                                                                                            ###
### Replication code for "L0 trend filtering" by Canhong Wen, Xueqin Wang and Aijun Zhang      ###
### This file is used to combine the RData and needed to be run before generating all the      ###
### figures in Appendix B.2.                                                                   ###
### Updated on 10 April 2023                                                                   ###
###                                                                                            ###
##################################################################################################



load("plot_results/tBlocks.RData")
ResultL0 <- ResultL0[1:8,]
ResultL1 <- ResultL1[1:8,]
ResultLc <- ResultLc[1:8,]
ResultNt <- ResultNt[1:8,]
ResultPt <- ResultPt[1:8,]
ResultSp <- ResultSp[1:8,]
ResultWs <- ResultWs[1:8,]
save(ResultL0, ResultL1, ResultLc, ResultNt, ResultSp, ResultPt, ResultWs, 
     file="pre_tBlocks.RData")
load("plot_results/tBlocks.RData")
ResultL0 <- ResultL0[8:15,]
ResultL1 <- ResultL1[8:15,]
ResultNt <- ResultNt[8:15,]
ResultPt <- ResultPt[8:15,]
ResultSp <- ResultSp[8:15,]
ResultWs <- ResultWs[8:15,]
save(ResultL0, ResultL1, ResultNt, ResultSp, ResultPt, ResultWs, 
     file="aft_tBlocks.RData")


load("plot_results/Blocks.RData")
ResultL0 <- ResultL0[1:8,]
ResultL1 <- ResultL1[1:8,]
ResultLc <- ResultLc[1:8,]
ResultNt <- ResultNt[1:8,]
ResultPt <- ResultPt[1:8,]
ResultSp <- ResultSp[1:8,]
ResultWs <- ResultWs[1:8,]
save(ResultL0, ResultL1, ResultLc, ResultNt, ResultSp, ResultPt, ResultWs, 
     file="pre_Blocks.RData")
load("plot_results/Blocks.RData")
ResultL0 <- ResultL0[8:15,]
ResultL1 <- ResultL1[8:15,]
ResultNt <- ResultNt[8:15,]
ResultPt <- ResultPt[8:15,]
ResultSp <- ResultSp[8:15,]
ResultWs <- ResultWs[8:15,]
save(ResultL0, ResultL1, ResultNt, ResultSp, ResultPt, ResultWs, 
     file="aft_Blocks.RData")


load("plot_results/tWave.RData")
ResultL0 <- ResultL0[1:8,]
ResultL1 <- ResultL1[1:8,]
ResultLc <- ResultLc[1:8,]
ResultNt <- ResultNt[1:8,]
ResultSp <- ResultSp[1:8,]
ResultCp <- ResultCp[1:8,]
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultSp,ResultCp,
     file="pre_tWave.RData")
load("plot_results/tWave.RData")
ResultL0 <- ResultL0[8:15,]
ResultL1 <- ResultL1[8:15,]
ResultNt <- ResultNt[8:15,]
ResultSp <- ResultSp[8:15,]
ResultCp <- ResultCp[8:15,]
save(ResultL0,ResultL1,ResultNt,ResultSp,ResultCp,
     file="aft_tWave.RData")


load("plot_results/Wave.RData")
ResultL0 <- ResultL0[1:8,]
ResultL1 <- ResultL1[1:8,]
ResultLc <- ResultLc[1:8,]
ResultNt <- ResultNt[1:8,]
ResultSp <- ResultSp[1:8,]
ResultCp <- ResultCp[1:8,]
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultSp,ResultCp,
     file="pre_Wave.RData")
load("plot_results/Wave.RData")
ResultL0 <- ResultL0[8:15,]
ResultL1 <- ResultL1[8:15,]
ResultNt <- ResultNt[8:15,]
ResultSp <- ResultSp[8:15,]
ResultCp <- ResultCp[8:15,]
save(ResultL0,ResultL1,ResultNt,ResultSp,ResultCp,
     file="aft_Wave.RData")


load("plot_results/tDoppler.RData")
ResultL0 <- ResultL0[1:8,]
ResultL1 <- ResultL1[1:8,]
ResultLc <- ResultLc[1:8,]
ResultNt <- ResultNt[1:8,]
ResultSp <- ResultSp[1:8,]
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultSp,
     file="pre_tDoppler.RData")
load("plot_results/tDoppler.RData")
ResultL0 <- ResultL0[8:15,]
ResultL1 <- ResultL1[8:15,]
ResultNt <- ResultNt[8:15,]
ResultSp <- ResultSp[8:15,]
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultSp,
     file="aft_tDoppler.RData")

load("plot_results/Doppler.RData")
ResultL0 <- ResultL0[1:8,]
ResultL1 <- ResultL1[1:8,]
ResultLc <- ResultLc[1:8,]
ResultNt <- ResultNt[1:8,]
ResultSp <- ResultSp[1:8,]
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultSp,
     file="pre_Doppler.RData")
load("plot_results/Doppler.RData")
ResultL0 <- ResultL0[8:15,]
ResultL1 <- ResultL1[8:15,]
ResultNt <- ResultNt[8:15,]
ResultSp <- ResultSp[8:15,]
save(ResultL0,ResultL1,ResultNt,ResultSp,
     file="aft_Doppler.RData")

load("plot_results/pre_Blocks.RData")
load("lcBlocks.RData")
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultPt,ResultSp,ResultWs,file="plot_results/pre_Blocks.RData")

load("plot_results/pre_Wave.RData")
load("lcWave.RData")
for(nn in seq(192,256,32)){
     ResultLc = rbind(ResultLc,c(nn,Inf,Inf,Inf,Inf,Inf,Inf,Inf))
}
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultSp,ResultCp,file="plot_results/pre_Wave.RData")

load("plot_results/pre_Doppler.RData")
load("lcDoppler.RData")
for(nn in seq(128,256,32)){
     ResultLc = rbind(ResultLc,c(nn,Inf,Inf,Inf))
}
save(ResultL0,ResultL1,ResultLc,ResultNt,ResultSp,file="plot_results/pre_Doppler.RData")
