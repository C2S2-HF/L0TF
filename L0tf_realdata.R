##################################################################################################
###                                                                                            ###
### Replication code for "L0 trend filtering" by Canhong Wen, Xueqin Wang and Aijun Zhang      ###
### This file contains the codes for real data application in PM 2.5 data                      ###
### Updated on 11 July 2022                                                                    ###
###                                                                                            ###
##################################################################################################


######################################################
### Load data                                      ###
######################################################
library(zoo)
DataX <- read.csv("air_hourly.csv", header = TRUE, na.strings = "N.A.")
DataX$Time = as.POSIXlt(paste(DataX$DATE,DataX$HOUR-1), format="%d/%m/%Y %H")
head(DataX)
y = zoo(DataX$FSP, order.by = DataX$Time)
y = na.omit(y)


######################################################
### Run L0-TF                                      ###
######################################################

library(AMIAS)
y1 = as.vector(y)
n = length(y1)
res =  samias(y1, D_type="tfq", q=1, rho = n^2.5,  kmax=110, adjust = FALSE, t.max = 5)
save(y, res, file="AirQualityResl0.RData")

######################################################
### Run L1-TF                                      ###
######################################################
library(genlasso)
y1 = as.vector(y)
n = length(y1)
resL1  <- trendfilter(y=y1, ord=1, maxsteps = 6000) 
save(y, resL1, file="AirQualityResl1.RData")


######################################################
### Plots                                          ###
######################################################
load("AirQualityResl0.RData")
load("AirQualityResl1.RData")

## Figure 19 ----------------------------------------
png("figs/HKAQ.png",  pointsize = 20, width=1000, height=700, res = 120)
op <- par(mfrow=c(1,1), mar=c(4, 4,3,2))
par(mfrow=c(1,1))
xtick <- as.Date(as.yearmon(2017 + 0.25+ seq(0,11)/12))  # The start date of each month
a <- which(strftime(index(y), format="%Y-%m-%d") %in% as.character(xtick))[c(1, seq(13, 267, 24))]
plot(as.vector(y), type="l", col="grey30", pch=20, cex= 0.8,  xaxt="n", 
     xlab="Time", ylab="PM2.5 Index", main="HK Central/Western PM2.5")   

axis(1,at=a,labels=c(format(xtick[1], "%Y%m"), format(xtick[2:9], "%m"), 
                     format(xtick[10], "%Y%m"), format(xtick[11:12], "%m")))
par(op)
dev.off()

png("figs/HKAQ_Path.png",  pointsize = 20, width=1000, height=700, res = 120)
op <- par(mfrow=c(1,1), mar=c(4,4,3,2))
plot(res, type="vpath", add.label = FALSE, k = NA)
title(main="Solution Path")
par(op)
dev.off()

beta <- res$alpha.all
png("figs/HKAQ_Fits.png",  pointsize = 20, width=1650, height=1200, res = 120)
op <- par(mfrow=c(2,2))
xtick <- as.Date(as.yearmon(2017 + 0.25+ seq(0,11)/12))  # The start date of each month
df <- c(10, 30, 60, 100)
for(i in 1:length(df)){
  par(mar=c(3,4,3,2))
  plot(as.numeric(y), type='p', col='grey', pch=19, cex= 0.6, xaxt="n", 
       main=paste0("df = ", df[i]), 
       ylab = "PM2.5 Index", xlab = "")
  a <- which(strftime(index(y), format="%Y-%m-%d") %in% as.character(xtick))[c(1, seq(13, 267, 24))]
  axis(1,at=a,labels=c(format(xtick[1], "%Y%m"), format(xtick[2:9], "%m"), 
                       format(xtick[10], "%Y%m"), format(xtick[11:12], "%m")))
  
  idx0 <- which(res$df==df[i])
  lines(beta[, idx0], col = "red", type='l', lty=1, lwd=2) ## plot the L0-TF estimate for a given df
  
  idx1 <- which(resL1$df==df[i])[1] # Choosing the one with the largest lambda value
  lines(resL1$beta[, idx1], col = "blue", type='l', lty=1, lwd=2) ## plot the L1-TF estimate for a given df
}
par(op)
dev.off()


## Figure 20 ------
df_zoom <- c(30, 60)
beta <- res$alpha.all
betaL1 <- resL1$beta

for(i in 1:2){
  ## print out the detected knots
  idx1 <- which(resL1$df==df_zoom[i])[1] # Choosing the one with the largest lambda value
  AL1 <- which(abs(diff(diff(betaL1[, idx1])))>1e-8)
  print(index(y)[AL1])
  
  idx0 <- which(res$df==df_zoom[i])
  AL0 <- which(abs(diff(diff(beta[, idx0])))>1e-8)
  print(index(y)[AL0])
}



day.zoom <- seq(from = as.Date("2017-08-23"), to = as.Date("2017-09-23"), by="day")
day.zoom <- as.character(day.zoom)
idx.zoom <- which(strftime(index(y), format="%Y-%m-%d") %in% day.zoom)

png("figs/HKAQ_Zoom.png", pointsize = 16, width=1150, height=450)
op <- par(mfrow=c(1,2))
for(i in 1:length(df_zoom)){
  par(mar=c(3,4,3,2))
  plot(as.numeric(y)[idx.zoom], type='p', col='grey', pch=19, cex= 0.6, xaxt="n", 
       main=paste0("df = ", df_zoom[i]), 
       ylab = "PM2.5 Index", xlab = "")
  
  a <- rep(0, length(day.zoom))
  for(j in 1:length(day.zoom)){
    a[j] <- which(strftime(index(y)[idx.zoom], format="%Y-%m-%d") %in% day.zoom[j])[1]
  }
  axis(1,at=a,labels=format(as.Date(day.zoom), "%m-%d"))
  
  idx0 <- which(res$df==df_zoom[i])
  lines(beta[idx.zoom, idx0], col = "red", type='l', lty=1, lwd=2) ## plot the L0-TF estimate for a given df
  
  idx1 <- which(resL1$df==df_zoom[i])[1] # Choosing the one with the largest lambda value
  lines(resL1$beta[idx.zoom, idx1], col = "blue", type='l', lty=1, lwd=2) ## plot the L1-TF estimate for a given df
  
  rect(6, -10, 16, 120, density = 30, border = "transparent", angle = -30) # HATO
  rect(96, -10, 104, 120, density = 30, border = "transparent", angle = -30) # PAKHAR
  
}
par(op)
dev.off()