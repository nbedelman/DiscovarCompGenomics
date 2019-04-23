library(dplyr)
library(ggplot2)

setwd("~/Dropbox")


FD="HeliconiusComparativeGenomics/FdComparison/fd_sim1Output.txt"

######## define functions #######
getSegment <- function (x) {
  segment <- strsplit(x,"_")[[1]][2]
  return(segment)
}


########## read in data ########
miyagi <- read.csv("HeliconiusComparativeGenomics/FdComparison/miyagi_sim1Output_3.csv",
                   header=F, col.names = c("segment","outgroup","branchLength","introProb","blank"))
miyagi <- select(miyagi, -blank)
miyagi$segment <- as.character(miyagi$segment)


fdData <- read.csv(FD)
fdData$segment <- as.character(unlist(lapply(as.character(fdData$scaffold),getSegment)))
fdData$segment <- as.character(as.numeric(fdData$segment)-1)
fdData <- subset(fdData, fdData$ABBA+fdData$BABA>10)

combined <- inner_join(miyagi,fdData)
########### make simple plots ##########
colors=brewer.pal(8,"Paired")

branchLengths <- ggplot() +
  geom_histogram(data=subset(combined, outgroup==3),aes(x=branchLength),fill=colors[2])+
  geom_histogram(data=subset(combined, outgroup==2),aes(x=branchLength),fill=colors[6])+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=branchLength),fill=colors[4])
branchLengths

introProbs <- ggplot() +
  geom_histogram(data=subset(combined, outgroup==3),aes(x=introProb),fill=colors[2])+
  geom_histogram(data=subset(combined, outgroup==2),aes(x=introProb),fill=colors[6])+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=introProb),fill=colors[4])
introProbs

BLvsIP <- ggplot(data=combined) +
  geom_point(aes(x=branchLength,y=introProb, col=as.factor(outgroup)))+
  scale_color_manual(values=c(colors[4],colors[6],colors[2]))
BLvsIP

BLvsFd <- ggplot(data=combined) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)))+
  scale_color_manual(values=c(colors[4],colors[6],colors[2]))
BLvsFd

BLvsIP <- ggplot(data=combined) +
  geom_point(aes(x=introProb,y=fd, col=as.factor(outgroup)))+
  scale_color_manual(values=c(colors[4],colors[6],colors[2]))
BLvsIP
