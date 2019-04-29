#Recombination by Miyagi-called introgression vs ILS with Fd statistic

library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(agricolae)
library(extrafont)
loadfonts(device="postscript")     
pdfFonts(Times =pdfFont("Times New Roman"))

setwd("~/Dropbox")

TREES="HeliconiusComparativeGenomics/FdComparison/eratoClade_HeraRef_5KB_allTrees.exist.found.csv"

RECRATE="HeliconiusComparativeGenomics/FdComparison/HeratoWindows.5KBAbutting.LD.bed"

FD="HeliconiusComparativeGenomics/FdComparison/fd_5kb_all.txt"

QUIBL="HeliconiusComparativeGenomics/FdComparison/eratoClade_HeraRef_5KB.allWindows.miyagi.csv"

out="HeliconiusComparativeGenomics/FdComparison/eratoClade_HeraRef_5KB.allWindows_"
########## Define Functions #############
getChrom <- function (x) {
  chrom <- strsplit(x,"_")[[1]][1]
  return(chrom)
}

getEnd <- function (x) {
  end <- strsplit(x,"_")[[1]][2]
  return(end)
}

getAlnLength <- function(x) {
  return(allData[which(allData$segment==x),]$alnLength)
}

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


########## read and format data #############

HeratoGenome <- read.delim("HeliconiusComparativeGenomics/datasets/Herato.transitions_fullChroms.bed", col.names=c("seqnames","start","end","chromosome"), header=F)
HeratoGenome$chr <-  unlist(strsplit(as.character(HeratoGenome$chromosome),"chr"))[c(FALSE,TRUE)]


trees <- read.csv(TREES, header=FALSE, col.names = c("segment","tree"))
trees$segment <- as.character(trees$segment)
trees$chrom <- as.character(unlist(lapply(trees$segment,getChrom)))
trees$end <- as.numeric(unlist(lapply(trees$segment,getEnd)))*5000
trees$start <- trees$end-4999

miyagi <- read.csv(QUIBL,header=F, col.names = c("segNumber","outgroup","branchLength","introProb","blank"))
miyagi <- select(miyagi, -blank)
taxa=c("H. erato","H. hecalesia", "H. telesiphe")

fd <- read.csv(FD)
colnames(fd)[1]="segment"
fd$segment=as.character(fd$segment)
fd$end <- as.numeric(unlist(lapply(fd$segment,getEnd)))*5000
fd$start <- fd$end-4999
fd$chrom <- as.character(unlist(lapply(fd$segment,getChrom)))

rateFile=read.delim(RECRATE, col.names=c("scaffold","start","end","segment","LD"))

#combined <- miyagi
combined <- cbind(trees,miyagi)
combined <- inner_join(combined,fd, by="segment")
combined <- inner_join(combined,rateFile, by="segment")
combined$LD <- as.numeric(as.character(combined$LD))
combined=combined[complete.cases(combined),]
combined$rateQuint <- as.factor(as.integer(cut(combined$LD, quantile(combined$LD, probs= seq(0,1,.2), na.rm=T), include.lowest=TRUE)))
combined <- subset(combined, ABBA+BABA>10)
combined$fd <- ifelse(combined$fd<0,0,combined$fd)
combinedSubset <- combined[sample(seq(1,nrow(combined)),size = 20000, replace=F),]

##categorize Trees 
#in order of Figure 2A
goodTrees <- c("tree7","tree9","tree3","tree2","tree35","tree5","tree8","tree34")

colors=brewer.pal(8,"Paired")
colorOrder <- c(colors[6],colors[2],colors[3],colors[4],colors[7],colors[5],colors[8],colors[1])



########### make simple plots ##########

branchLengths <- ggplot() +
  geom_histogram(data=subset(combined, outgroup==3),aes(x=branchLength),fill=colors[2], bins=50)+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=branchLength),fill=colors[6], bins=50, alpha=0.5)+
  geom_histogram(data=subset(combined, outgroup==2),aes(x=branchLength),fill=colors[4], bins=50, alpha=0.5)+
  scale_x_continuous(limits=c(0,0.07))
branchLengths
ggsave(branchLengths,file=paste0(out,"branchLengthsHist.pdf"),height=10,width=10)

branchLengthBox <- ggplot(data=combinedSubset,aes(x=as.factor(outgroup),y=branchLength,col=as.factor(outgroup))) +
  geom_jitter(alpha=.1)+
  geom_boxplot(outlier.shape = NA,fill="transparent", color="black")+
  scale_color_manual(values=c(colors[6],colors[4],colors[2]), labels=taxa,name="outgroup")+
  scale_x_discrete(labels=taxa) +
  labs(x="outgroup")+
  scale_y_continuous(limits=c(0,0.07))+
  theme(legend.position = "none")
branchLengthBox
ggsave(branchLengthBox,file=paste0(out,"branchLengthsBox.pdf"),height=10,width=10)


introProbs <- ggplot() +
  geom_histogram(data=subset(combined, outgroup==3),aes(x=introProb),fill=colors[2],alpha=.7)+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=introProb),fill=colors[6],alpha=.7)+
  geom_histogram(data=subset(combined, outgroup==2),aes(x=introProb),fill=colors[4],alpha=.7)+
  labs(x="Introgression Probability",y="Number of Windows")
introProbs
ggsave(introProbs,file=paste0(out,"introgressionProbHist.pdf"),height=10,width=10)


BLvsIP <- ggplot(data=combined[sample(seq(1,nrow(combined)),size = 5000, replace=F),]) +
  geom_point(aes(x=branchLength,y=introProb, col=as.factor(outgroup)), alpha=.2)+
  scale_color_manual(values=c(colors[6],colors[4],colors[2]), labels=taxa)+
  scale_x_continuous(limits=c(0,0.07))+
  theme(legend.position = "none")#+
  #geom_vline(xintercept =.01)
BLvsIP
ggsave(BLvsIP,file=paste0(out,"BranchLengthvsIntroProb.pdf"),height=10,width=10)


BLvsFd <- ggplot(data=combined) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.7)+
  scale_color_manual(values=c(colors[6],colors[4],colors[2]))+
  scale_x_continuous(limits=c(0,0.1))
BLvsFd
ggsave(BLvsFd,file=paste0(out,"BranchLengthvsFd.pdf"),height=10,width=10)


BLvsFd_intro <- ggplot(data=subset(combinedSubset,outgroup==1)) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.2)+
  scale_color_manual(values=c(colors[6]))+
  scale_x_continuous(limits=c(0,0.07))+
  theme(legend.position = "none") 
BLvsFd_intro
ggsave(BLvsFd_intro,file=paste0(out,"BranchLengthvsFd_intro.pdf"),height=10,width=10)


BLvsFd_species <- ggplot(data=subset(combined,outgroup==3)) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.1)+
  scale_color_manual(values=c(colors[2]))+
  scale_x_continuous(limits=c(0,0.1))
BLvsFd_species
ggsave(BLvsFd_species,file=paste0(out,"BranchLengthvsFd_species.pdf"),height=10,width=10)

BLvsFd_ILS <- ggplot(data=subset(combined,outgroup==2)) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.1)+
  scale_color_manual(values=c(colors[4]))+
  scale_x_continuous(limits=c(0,0.1))
BLvsFd_ILS
ggsave(BLvsFd_ILS,file=paste0(out,"BranchLengthvsFd_ILS.pdf"),height=10,width=10)

IPvsFd_ILS <- ggplot(data=subset(combined,outgroup==2)) +
  geom_point(aes(x=introProb,y=fd, col=as.factor(outgroup)), alpha=.1)+
  scale_color_manual(values=c(colors[4]))+
  geom_smooth(aes(x=introProb,y=fd), col="black")
IPvsFd_ILS
ggsave(IPvsFd_ILS,file=paste0(out,"IntroProbvsFd_ILS.pdf"),height=10,width=10)


IPvsFd_intro <- ggplot(data=subset(combinedSubset,outgroup==1)) +
  geom_point(aes(x=introProb,y=fd, col=as.factor(outgroup)), alpha=.2)+
  scale_color_manual(values=c(colors[6]))+
  geom_smooth(aes(x=introProb,y=fd), col="black")+
  #scale_y_continuous(limits=c(0,0.65))+
  theme(legend.position="none")
IPvsFd_intro
ggsave(IPvsFd_intro,file=paste0(out,"IntroProbvsFd_intro.pdf"),height=10,width=10)

Fd_maxIP <- ggplot(data=subset(combined,outgroup==1 & introProb >= max(subset(combined,outgroup==1)$introProb)))+
  geom_density(aes(x=fd),fill="red")+
  geom_density(aes(x=fd),fill="blue", alpha=.75)
Fd_maxIP

FdVsFd <- ggplot(data=combined)+
  geom_density(aes(x=fdM),fill="red")+
  geom_density(aes(x=fd),fill="blue", alpha=.75)
FdVsFd


##### Plot mean Fd by chromosome size ########

meanfd <- c()
meanfd_intro <- c()
meanfd_withNeg <- c()
meanfd_abs <- c()
meanfd_intro_abs <- c()
meanIP <- c()
meanBL <- c()
for (chr in HeratoGenome$chromosome){
  meanfd <- c(meanfd,mean(combined[which(combined$chrom.x == chr & combined$fd>=0),]$fd))
  meanfd_withNeg <- c(meanfd_withNeg,mean(combined[which(combined$chrom.x == chr),]$fd))
  meanfd_abs <- c(meanfd_abs,mean(abs(combined[which(combined$chrom.x == chr),]$fd)))
  meanfd_intro <- c(meanfd_intro,mean(subset(combined,outgroup==1 & fd>=0)[which(subset(combined,outgroup==1 & fd>=0)$chrom.x == chr),]$fd))
  meanfd_intro_abs<- c(meanfd_intro_abs,mean(abs(subset(combined,outgroup==1)[which(subset(combined,outgroup==1)$chrom.x == chr),]$fd)))
  meanIP <- c(meanIP,mean(subset(combined,outgroup==1)[which(subset(combined,outgroup==1)$chrom.x == chr),]$introProb))
  meanBL <- c(meanBL,mean(subset(combined,outgroup==1)[which(subset(combined,outgroup==1)$chrom.x == chr),]$branchLength))
}

HeratoGenome$meanfd <- meanfd
HeratoGenome$meanfd_intro <- meanfd_intro
HeratoGenome$meanfd_intro_abs <- meanfd_intro_abs
HeratoGenome$meanfd_abs <- meanfd_abs
HeratoGenome$meanIP <- meanIP
HeratoGenome$meanBL <- meanBL
HeratoGenome$lengthMB <- HeratoGenome$end/1000000
length_rho <- cor.test(~meanfd+lengthMB, data=subset(HeratoGenome, chr !=21), method="spearman")
length_lm <- lm(meanfd_abs~lengthMB, data=subset(HeratoGenome, chr !=21))
colors=brewer.pal(8,"Paired")
fdByLength <- ggplot(data=HeratoGenome) +
  geom_point(aes(y=meanfd,x=lengthMB),col=colors[6]) +
  geom_text_repel(aes(label=chromosome,y=rep(max(meanfd+.005),21),x=lengthMB),direction="x",nudge_y=.01)+
  geom_abline(slope=length_lm$coefficients[2],intercept=length_lm$coefficients[1],col=colors[6])+
  # #ylim(0,0.70)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Length (MB)",y=expression(italic("f"['d'])))
fdByLength
ggsave(fdByLength,file=paste0(out,"fdvsChromLength.pdf"),height=10,width=10)


length_fd_intro <- lm(meanfd_intro~lengthMB, data=subset(HeratoGenome, chr !=21))
colors=brewer.pal(8,"Paired")
fd_introByLength <- ggplot(data=HeratoGenome) +
  geom_point(aes(y=meanfd_intro,x=lengthMB),col=colors[6]) +
  geom_text_repel(aes(label=chromosome,y=rep(max(meanfd_intro+.005),21),x=lengthMB),direction="x",nudge_y=.01)+
  geom_abline(slope=length_fd_intro$coefficients[2],intercept=length_fd_intro$coefficients[1],col=colors[6])+
  # #ylim(0,0.70)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Length (MB)",y=expression(italic("f"['d'])))
fd_introByLength
ggsave(fd_introByLength,file=paste0(out,"fdvsChromLength_introTop.pdf"),height=10,width=10)

IPLength_lm <- lm(meanIP~lengthMB, data=subset(HeratoGenome, chr !=21))
IPbyLength <- ggplot(data=HeratoGenome) +
  geom_point(aes(y=meanIP,x=lengthMB),col=colors[6]) +
  geom_text_repel(aes(label=chromosome,y=rep(max(meanIP+.002),21),x=lengthMB),direction="x",nudge_y=.002)+
  geom_abline(slope=IPLength_lm$coefficients[2],intercept=IPLength_lm$coefficients[1],col=colors[6])+
  ylim(0.94,0.958)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Length (MB)",y="Introgression Probability")
IPbyLength
ggsave(IPbyLength,file=paste0(out,"IntroProbvsChromLength.pdf"),height=10,width=10)


BLLength_lm <- lm(meanBL~lengthMB, data=subset(HeratoGenome, chr !=21))
BLbyLength <- ggplot(data=HeratoGenome) +
  geom_point(aes(y=meanBL,x=lengthMB),col=colors[6]) +
  geom_text_repel(aes(label=chromosome,y=rep(max(meanBL+.001),21),x=lengthMB),direction="x",nudge_y=.002)+
  geom_abline(slope=BLLength_lm$coefficients[2],intercept=BLLength_lm$coefficients[1],col=colors[6])+
  # #ylim(0,0.70)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Length (MB)",y="Branch Length")
BLbyLength
ggsave(BLbyLength,file=paste0(out,"BranchLengthvsChromLength.pdf"),height=10,width=10)

##### plot Fd and intro prob by local recombination rate #######

#unbinned
fdModel <- lm(fd~LD, data=combined)
lmIntercept=fdModel$coefficients[1]
lmSlope=fdModel$coefficients[2]
confFrame <- data.frame(LD=combined$LD)
confInterval <- as.data.frame(predict(fdModel, confFrame, interval="confidence"))
fdVsRate <- ggplot(data=combined)+
  geom_point(aes(x=LD,y=fd), alpha=0.05)+
  geom_line(aes(x=confFrame$LD, y=confInterval$lwr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$upr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$fit), lty=1, col="orange")+
  #geom_smooth(aes(x=LD,y=fd))+
  labs(x="recombination rate",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
fdVsRate
ggsave(fdVsRate,file=paste0(out,"FdvsRecRate_all.pdf"),height=10,width=10)

###use only the introgression topology
onlyIntro <- subset(combined,outgroup==1)
onlyIntro$rateQuint <- as.factor(as.integer(cut(onlyIntro$LD, quantile(onlyIntro$LD, probs= seq(0,1,.2), na.rm=T), include.lowest=TRUE)))

fdModel <- lm(fd~LD, data=onlyIntro)
lmIntercept=fdModel$coefficients[1]
lmSlope=fdModel$coefficients[2]
confFrame <- data.frame(LD=onlyIntro$LD)
confInterval <- as.data.frame(predict(fdModel, confFrame, interval="confidence"))
fdVsRate_intro <- ggplot(data=onlyIntro)+
  geom_point(aes(x=LD,y=fd), alpha=0.05)+
  geom_line(aes(x=confFrame$LD, y=confInterval$lwr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$upr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$fit), lty=1, col="orange")+
  #geom_smooth(aes(x=LD,y=fd))+
  labs(x="recombination rate",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
fdVsRate_intro
ggsave(fdVsRate_intro,file=paste0(out,"FdvsRecRate_introTop.pdf"),height=10,width=10)


IPModel <- lm(introProb~LD, data=onlyIntro)
lmIntercept=IPModel$coefficients[1]
lmSlope=IPModel$coefficients[2]
confFrame <- data.frame(LD=onlyIntro$LD)
confInterval <- as.data.frame(predict(IPModel, confFrame, interval="confidence"))
IPVsRate <- ggplot(data=onlyIntro)+
  geom_point(aes(x=LD,y=introProb), alpha=0.05)+
  geom_line(aes(x=confFrame$LD, y=confInterval$lwr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$upr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$fit), lty=1, col="orange")+
  #geom_smooth(aes(x=LD,y=introProb))+
  labs(x="recombination rate",y="Introgression Probability")+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
IPVsRate
ggsave(IPVsRate,file=paste0(out,"IntroProbvsRecRate_introTop.pdf"),height=10,width=10)



BLModel <- lm(branchLength~LD, data=onlyIntro)
lmIntercept=BLModel$coefficients[1]
lmSlope=BLModel$coefficients[2]
confFrame <- data.frame(LD=onlyIntro$LD)
confInterval <- as.data.frame(predict(BLModel, confFrame, interval="confidence"))
BLVsRate <- ggplot(data=onlyIntro)+
  geom_point(aes(x=LD,y=branchLength), alpha=0.05)+
  geom_line(aes(x=confFrame$LD, y=confInterval$lwr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$upr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$fit), lty=1, col="orange")+
  #geom_smooth(aes(x=LD,y=branchLength))+
  labs(x="recombination rate",y="Internal Branch Length")+
  scale_y_continuous(limits=c(0,0.06))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
BLVsRate
ggsave(BLVsRate,file=paste0(out,"BranchLengthvsRecRate_introTop.pdf"),height=10,width=10)


##### separate into bins ########
##all data
rateKW <- kruskal(combined$fd, combined$rateQuint, group=TRUE, p.adj="bonferroni")
ratedf <- data.frame(tree=rownames(rateKW$groups), rateQuint=rateKW$groups$groups)

FdByRate <-   ggplot(data=combined)+
  geom_jitter(aes(y=fd,x=rateQuint), alpha=.05)+
  geom_boxplot(aes(y=fd, x=rateQuint),outlier.shape = NA,col=colors[6], lwd=1, fill="transparent")+
  geom_text(data=ratedf, aes(x=tree, y=0.8, label=rateQuint), size=8)+
  labs(x="recombination rate quintile",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
FdByRate
ggsave(FdByRate,file=paste0(out,"BranchLengthvsRecRate_quint_all.pdf"),height=10,width=10)

###introgression topology
fdRate_introKW <- kruskal(onlyIntro$fd, onlyIntro$rateQuint, group=TRUE, p.adj="bonferroni")
fdRate_introdf <- data.frame(tree=rownames(fdRate_introKW$groups), rateQuint=fdRate_introKW$groups$groups)

FdByRate_intro <-   ggplot(data=onlyIntro)+
  geom_jitter(aes(y=fd,x=rateQuint), alpha=.05)+
  geom_boxplot(aes(y=fd, x=rateQuint),outlier.shape = NA,col=colors[6], lwd=1, fill="transparent")+
  geom_text(data=fdRate_introdf, aes(x=tree, y=0.8, label=rateQuint), size=8)+
  labs(x="recombination rate quintile",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
FdByRate_intro
ggsave(FdByRate_intro,file=paste0(out,"FdvsRecRate_quint_introTop.pdf"),height=10,width=10)

IPRate_introKW <- kruskal(onlyIntro$introProb, onlyIntro$rateQuint, group=TRUE, p.adj="bonferroni")
IPRate_introdf <- data.frame(tree=rownames(IPRate_introKW$groups), rateQuint=IPRate_introKW$groups$groups)

IPByRate <-   ggplot(data=onlyIntro)+
  geom_jitter(aes(y=introProb,x=rateQuint), alpha=.05)+
  geom_boxplot(aes(y=introProb, x=rateQuint),outlier.shape = NA,col=colors[6], lwd=1, fill="transparent")+
  geom_text(data=IPRate_introdf, aes(x=tree, y=1, label=rateQuint), size=8)+
  labs(x="recombination rate quintile",y="Introgression probability")+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
IPByRate
ggsave(IPByRate,file=paste0(out,"IntroProbvsRecRate_quint_introTop.pdf"),height=10,width=10)

BLRate_introKW <- kruskal(onlyIntro$branchLength, onlyIntro$rateQuint, group=TRUE, p.adj="bonferroni")
BLRate_introdf <- data.frame(tree=rownames(BLRate_introKW$groups), rateQuint=BLRate_introKW$groups$groups)

BLByRate <-   ggplot(data=onlyIntro)+
  geom_jitter(aes(y=branchLength,x=rateQuint), alpha=.05)+
  geom_boxplot(aes(y=branchLength, x=rateQuint),outlier.shape = NA,col=colors[6], lwd=1, fill="transparent")+
  geom_text(data=BLRate_introdf, aes(x=tree, y=0.045, label=rateQuint), size=8)+
  labs(x="recombination rate quintile",y="Branch Length")+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))+
  scale_y_continuous(limits=c(0,0.05))
BLByRate
ggsave(IPByRate,file=paste0(out,"BranchLengthvsRecRate_quint_introTop.pdf"),height=10,width=10)



######### Fd by tree ########
treeSubsetFrame <- data.frame()
for (treeType in goodTrees){
  treeSub <- combined[sample(which(combined$tree==treeType),750, replace = F),]
  treeSubsetFrame <- rbind(treeSubsetFrame,treeSub )
}
treeSubsetFrame$tree <- factor(treeSubsetFrame$tree, levels=goodTrees)

treeKW <- kruskal(treeSubsetFrame$fd, treeSubsetFrame$tree, group=TRUE, p.adj="bonferroni")
treedf <- data.frame(tree=rownames(treeKW$groups), treeGroup=treeKW$groups$groups)
treedf$tree <- factor(treedf$tree, levels=goodTrees)

FdByTree <-   ggplot(data=subset(treeSubsetFrame, outgroup==1))+
  geom_boxplot(aes(y=fd, x=tree),outlier.shape = NA)+
  geom_jitter(aes(y=fd,x=tree), alpha=.05)+
  geom_text(data=treedf, aes(x=tree, y=0.6, label=treeGroup), size=8)+
  labs(x="",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))+
  scale_x_discrete(limits=goodTrees,labels=c("Tree 1","Tree 2","Tree 3","Tree 4","Tree 5","Tree 6","Tree 7","Tree 8"))
FdByTree

####branch length by tree
BLByTree <-   ggplot(data=subset(treeSubsetFrame, outgroup==1))+
  geom_boxplot(aes(y=branchLength, x=tree),outlier.shape = NA)+
  geom_jitter(aes(y=branchLength,x=tree), alpha=.05)+
  scale_x_discrete(labels=c("Tree 1","Tree 2","Tree 3","Tree 4","Tree 5","Tree 6","Tree 7","Tree 8"),limits=goodTrees)+
  #geom_text(data=treedf, aes(x=tree, y=0.6, label=treeGroup), size=8)+
  labs(x="",y="Internal Branch Length")+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
BLByTree

####Intro Prob by tree
IPByTree <-   ggplot(data=subset(treeSubsetFrame, outgroup==1))+
  geom_boxplot(aes(y=introProb, x=tree),outlier.shape = NA)+
  geom_jitter(aes(y=introProb,x=tree), alpha=.05)+
  scale_x_discrete(labels=c("Tree 1","Tree 2","Tree 3","Tree 4","Tree 5","Tree 6","Tree 7","Tree 8"),limits=goodTrees)+
  #geom_text(data=treedf, aes(x=tree, y=0.6, label=treeGroup), size=8)+
  labs(x="",y="Internal Branch Length")+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
IPByTree




####### D by tree #########
DByTree <-   ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=D, x=tree),outlier.shape = NA)+
  geom_jitter(aes(y=D,x=tree), alpha=.05)+
  scale_x_discrete(labels=c("Tree 1","Tree 2","Tree 3","Tree 4","Tree 5","Tree 6","Tree 7","Tree 8"))+
  #geom_text(data=treedf, aes(x=tree, y=0.6, label=treeGroup), size=8)+
  labs(x="",y=expression(italic("D")))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
DByTree


####### make figure #######

ng <- nullGrob()
g <- arrangeGrob(FdByRate, FdByTree, 
                 FdByCall, ng)
ggsave(filename = "fd_recomb_miyagi5KBevery50_kwtest.pdf",g, width = 10, height=10)
