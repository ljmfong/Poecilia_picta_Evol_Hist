#cov_fc_mov_av

rm(list=ls())
ls() 

library(ggplot2)
library(reshape)
library(gridExtra)
library(dplyr)
library(ggthemes)
library(extrafont)
loadfonts()
library(RcppRoll)
library(zoo)

#define a moving average function
movingaverage <- function (x, window) {
  ma <- roll_mean(x, window, fill = NA)
}
windowsizecov = 50

#coverage data
#get data
datacov <- read.csv(file.choose(), stringsAsFactors=FALSE,  header=T) #This should also have positional info
autocov <- read.csv(file.choose(), stringsAsFactors = FALSE, header = T)
#This will be your MF_foldchange.csv file 
dim(datacov)
# 6259   8
dim(autocov)
# 300   10
names(datacov)
# "ï..scaffold" "binstart"    "binend"      "MFlogavg"    "Mavg"        "Mlogavg"     "Favg"        "Flogavg"
names(autocov)

colnames(datacov)[1] <- "scaffold" 
colnames(autocov)[1] <- "scaffold"

#Double check to make sure it's been renamed properly
names(datacov)

mean(datacov$MFlogavg)
#-0.7123898

mean(autocov$Mflogavg)
#-0.02755158

#sort data - in this case, the chromosome is already sorted into the bins of 50kB
autocovsorted <- autocov[order(autocov$scaffold,autocov$start),]
dim(autocovsorted)
# 300   10

#calculate confidence interval for moving average

MFautopermute <- replicate(1000,mean(sample(autocovsorted$Mflogavg,windowsizecov,replace = FALSE)))
MFautoI25cov <- quantile(MFautopermute, c(.025, .5, .975))[[1]]
MFautoI25cov
#-0.1767492
MFautoI975cov <- quantile(MFautopermute, c(.025, .5, .975))[[3]]
MFautoI975cov
#0.1338285

MFpermutescov <- replicate(1000,mean(sample(datacov$MFlogavg,windowsizecov,replace = FALSE)))
MFCI25cov <- quantile(MFpermutescov, c(.025, .5, .975))[[1]]
MFCI25cov
# -0.8239097
MFCI975cov <- quantile(MFpermutescov, c(.025, .5, .975))[[3]]
MFCI975cov
# -0.5757917

#This was when I was only using PAR
datacovAUTOSOME <- data.frame(datacov$binstart[6000:6259],datacov$MFlogavg[6000:6259])
colnames(datacovAUTOSOME)[1] <- "PAR"
colnames(datacovAUTOSOME)[2] <- "MFlogavg"


MFpermutescovPAR <- replicate(1000,mean(sample(datacovAUTOSOME$MFlogavg,windowsizecov,replace = FALSE)))
MFCI25covPAR <- quantile(MFpermutescovPAR, c(.025, .5, .975))[[1]]
MFCI25covPAR
# -0.1872198
MFCI975covPAR <- quantile(MFpermutescovPAR, c(.025, .5, .975))[[3]]
MFCI975covPAR
# 0.04079718
#

#plot chr 12 coverage data

median(datacovAUTOSOME$MFlogavg)
# -0.1685105

#plot chr 12 snp data

datacovLG12 <- datacov[order(datacov$binstart,as.numeric(levels(datacov$binstart))[datacov$binstart]),]

dim(datacovLG12)

smoothlinefccov = movingaverage(datacov$MFlogavg,windowsizecov)
smoothlinefemcov = movingaverage(datacov$Flogavg,windowsizecov)
smoothlinemalcov = movingaverage(datacov$Mlogavg,windowsizecov)


paneldatacov <- as.data.frame(smoothlinefccov)
paneldatacov$startmb <- datacov$binstart
paneldatacov$females <- smoothlinefemcov
paneldatacov$males <- smoothlinemalcov
names(paneldatacov)
# "smoothlinefccov" "startmb"         "females"         "males"
paneldatacov <- na.omit(paneldatacov)

# plots (you might have to change the axes to view your plot better - depends on the values plotted)
require(grid)

plotfccov <- ggplot(datacov, aes(x= binstart/1000000, y= MFlogavg)) +
  geom_rect(xmax=60,xmin=-10,ymax=MFautoI975cov,ymin=MFautoI25cov,fill="grey", alpha=0.08)+
  geom_point(colour="grey50",fill="grey50", alpha=0.5, cex=0.4) +
    geom_line(data = paneldatacov, aes(x=  paneldatacov$startmb/1000000, y= paneldatacov$smoothlinefccov), colour = "black", size=0.6) +
    coord_cartesian(ylim=c(-2,2)) +
    coord_cartesian(xlim=c(0,31.5)) +
    theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
    theme(
          text=element_text(size=8),
          plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
          axis.line.y = element_line(color="black", size = 0.3),
          axis.line.x = element_line(color="black", size = 0.3),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8)
          ) +
      ylim(-2,2) + 
      scale_x_continuous(breaks=seq(0,35,5)) +
      xlab('Start position (Mb)') +
      ylab(expression('M:F log'[2]*' coverage')) +
  annotate(geom = "rect", xmin=29.995,xmax=31.5,ymin=-2, ymax=2, 
            fill = "darkorange3", alpha=0.2)

plotfccov



plotsexcov <- ggplot(datacov, aes(x= binstart/1000000, y= Flogavg)) + 
  geom_point(colour="#e31a1c",fill="#e31a1c", alpha=0.5, cex=0.4) +
  geom_point(data= datacov, aes(x= binstart/1000000, y= Mlogavg),colour="#1f78b4",fill="#1f78b4", alpha=0.5, cex=0.4) +
  geom_line(data = paneldatacov, aes(x= paneldatacov$startmb/1000000, y=paneldatacov$females), colour="#e31a1c",, size=0.6) +
  geom_line(data = paneldatacov, aes(x= paneldatacov$startmb/1000000, y=paneldatacov$males), colour="#1f78b4",, size=0.6) +
  coord_cartesian(ylim=c(5,7)) +
  coord_cartesian(xlim=c(0,31.5)) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
  theme(
    text=element_text(size=8),
    plot.margin = unit(c(1.5,1.5,1.5,1.5),"lines"),
    axis.line.y = element_line(color="black", size = 0.3),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.text.x = element_text(size=8),
    axis.text.y = element_text(size=8)
  ) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks=seq(0,35,5)) +
  xlab('Start position (Mb)') +
  ylab(expression('Log'[2]*' coverage'))

plotsexcov

require(cowplot)
plot_grid(plotfccov, plotsexcov, labels = c('A', 'B'),ncol = 1,nrow=2, align="v")
