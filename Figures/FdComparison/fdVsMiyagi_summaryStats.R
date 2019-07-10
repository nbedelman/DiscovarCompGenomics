### compare general outputs from Fd and QuIBL

library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(agricolae)

setwd("~/Dropbox/HeliconiusComparativeGenomics/FdComparison/")

#####define datafiles 
QUIBL2DIST="miyagi_outputSummary_all2Dist.csv"
QUIBL1DIST="miyagi_outputSummary_all1Dist.csv"
FdSummary="eratoClade.helSingleCopy.Fd_withNeg.csv"

out="miyagi_outputSummary_allData_"


#####read in data 
quibl2Dist=read.csv(QUIBL2DIST,header=F,col.names=c("triplet","P1","C1","C2","prop1","prop2","numTrees","BIC","maxL"))
quibl1Dist=read.csv(QUIBL1DIST,header=F,col.names=c("triplet","P1","C1","C2","prop1","prop2","numTrees","BIC","maxL"))


quibl <- cbind(quibl2Dist,quibl1Dist$BIC,quibl1Dist$maxL)
colnames(quibl)[10] <- "BIC1"
colnames(quibl)[11] <- "maxL1"
quibl$BICdiff=quibl$BIC-quibl$BIC1
quibl$likeDiff=quibl$maxL-quibl$maxL1
quibl$totalIntroProp <- quibl$prop2*(quibl$numTrees/5600)
quibl$isSig <- quibl$BICdiff < (-15)
quibl$MLsig <- quibl$likeDiff> (-10)

fd=read.csv(FdSummary)
triplets <- c()
for(line in seq(1,nrow(fd))){
  seqList=c(as.character(fd[line,]$P1),as.character(fd[line,]$P2),as.character(fd[line,]$P3))
  seqList <- sort(seqList)
  triplets <- c(triplets,paste(seqList[1],seqList[2],seqList[3],sep="_"))
}
fd$triplet <- triplets




outGroup="HmelRef"
if (is.na(outGroup)){
  quibl$hasOverallOut <- F
} else {quibl$hasOverallOut <- grepl(outGroup,quibl$triplet)}
quibl <- subset(quibl,hasOverallOut==F)

combined <- inner_join(quibl,fd,by=c("triplet","P1"))

geneTreeIndices=c(1,3,7,9,13,15,19,21,25,27,31,33,35,39,41,45,47,54,56,57,59,63,65,69,71,75,77,81,83,90,92,93,95,99,101,108,110,114,116)
fd_geneTree <- fd[geneTreeIndices,]
combined_geneTree <- inner_join(quibl,fd_geneTree,by=c("triplet","P1"))

colors=brewer.pal(8,"Paired")

##### total introgression proportions #####
#here, we're calculating the total number of discordant triplets that arose via introgression. Assuming that concordant triplets are most common.
#EtalOut 3330/2
#HQ 1800
#allData 2790
discordSet=subset(quibl, numTrees<(2260))
allDiscordTrees=sum(discordSet$numTrees)
pctDiscordTrees=allDiscordTrees/sum(quibl$numTrees)
discordIntro <- sum(discordSet$numTrees*discordSet$prop2*as.numeric(discordSet$isSig))
#discordIntro <- sum(discordSet$numTrees*discordSet$prop2*as.numeric(discordSet$MLsig))
#discordIntro <- sum(discordSet$numTrees*discordSet$prop2)
propDiscordIntro <- discordIntro/allDiscordTrees
propOverallIntro <- discordIntro/sum(quibl$numTrees)

length(which(discordSet$isSig==T))
length(unique(subset(discordSet, isSig==T)$triplet))

#### overview graphs 
CvalHist <- ggplot(data=combined,aes(x=C2))+
  geom_histogram(bins=20)+
  labs(y="Number of topologies",x="introgression C")
  #geom_vline(xintercept = 0.5)
CvalHist
ggsave(CvalHist,file=paste0(out,"CVal_hist.pdf"),height=10,width=10)

introPropHist <- ggplot(data=combined,aes(x=prop2))+
  geom_histogram(bins=20)+
  labs(y="Number of topologies",x="Topology introgression proportion")
introPropHist
ggsave(introPropHist,file=paste0(out,"topIntroProp_hist.pdf"),height=10,width=10)

CvalVsIntroProp <- ggplot(data=combined)+
  geom_point(aes(x=C2,y=prop2,size=numTrees/5600, col=numTrees>2790, pch=isSig), alpha=.7)+
  labs(y="Topology introgression proportion",x="Introgression C")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minority","majority"),name="")+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")+
  theme(legend.position="none")
CvalVsIntroProp
ggsave(CvalVsIntroProp,file=paste0(out,"CValVsTopIntroProp.pdf"),height=10,width=10)

CvalVsTotalIntroProp <- ggplot(data=combined)+
  geom_point(aes(x=C2,y=totalIntroProp,size=numTrees/5600, col=numTrees>2790, pch=isSig), alpha=.7)+
  labs(y="Total introgression proportion",x="Introgression C")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minority","majority"),name="")+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")+
  theme(legend.position="none")
CvalVsTotalIntroProp
ggsave(CvalVsTotalIntroProp,file=paste0(out,"CValVsTotalIntroProp.pdf"),height=10,width=10)

FdVsTotalIntroProp <- ggplot(data=combined_geneTree)+
  geom_point(aes(x=totalIntroProp,y=meanFd,size=numTrees/5600, col=numTrees>2790, pch=isSig), alpha=.7)+
  labs(x="Total introgression proportion",y="meanFd")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minority","majority"),name="")+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")+
  theme(legend.position="none")+
  scale_x_continuous(limits=c(0,0.5))+
  scale_y_continuous(limits=c(0,0.5))+
  geom_abline(slope=1,intercept=0)
FdVsTotalIntroProp
ggsave(FdVsTotalIntroProp,file=paste0(out,"FdVsTotalIntroProp_geneTrees.pdf"),height=10,width=10)

DVsTotalIntroProp <- ggplot(data=combined_geneTree)+
  geom_point(aes(x=totalIntroProp,y=abs(meanD),size=numTrees/5600, col=numTrees>2790, pch=isSig), alpha=.7)+
  labs(x="Total introgression proportion",y="abs(mean D)")+
  scale_size_continuous(name = "Proportion of topologies")+
  scale_color_manual(values=c(colors[8],colors[2]),labels=c("minority","majority"),name="")+
  scale_x_continuous(limits=c(0,0.5))+
  scale_y_continuous(limits=c(0,0.5))+
  geom_abline(slope=1,intercept=0)+
  scale_shape_manual(values=c(18,19),labels=c("nonsignificant","significant"),name="")+
  theme(legend.position="none")
DVsTotalIntroProp
ggsave(DVsTotalIntroProp,file=paste0(out,"DVsTotalIntroProp_geneTrees.pdf"),height=10,width=10)

