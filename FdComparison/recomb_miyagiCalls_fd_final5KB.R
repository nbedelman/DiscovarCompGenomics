#Recombination by Miyagi-called introgression vs ILS with Fd statistic

library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
#library(agricolae)

setwd("~/Dropbox")

TREES="HeliconiusComparativeGenomics/FdComparison/eratoClade_HeraRef_5KB_highQual.found"

RECRATE="HeliconiusComparativeGenomics/FdComparison/HeratoWindows.5KBAbutting.LD.bed"

FD="HeliconiusComparativeGenomics/FdComparison/eratoClade_HeraRef_5KB_pRecombLe.05.miyagi.csv"

QUIBL="HeliconiusComparativeGenomics/FdComparison/eratoClade_HeraRef_5KB_pRecombLe.05.miyagi.csv"
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

fd <- read.csv(FD)
colnames(fd)[1]="segment"
fd$segment=as.character(fd$segment)
fd$end <- as.numeric(unlist(lapply(fd$segment,getEnd)))*5000
fd$start <- fd$end-4999
fd$chrom <- as.character(unlist(lapply(fd$segment,getChrom)))

rateFile=read.delim(RECRATE, col.names=c("scaffold","start","end","segment","LD"))


combined <- cbind(trees,miyagi)
combined <- inner_join(combined,fd, by="segment")
combined <- inner_join(combined,rateFile, by="segment")
combined$LD <- as.numeric(as.character(combined$LD))
combined=combined[complete.cases(combined),]
combined$rateQuint <- as.factor(as.integer(cut(combined$LD, quantile(combined$LD, probs= seq(0,1,.2), na.rm=T), include.lowest=TRUE)))

########### make simple plots ##########
colors=brewer.pal(8,"Paired")

branchLengths <- ggplot() +
  #geom_histogram(data=subset(combined, outgroup==3),aes(x=branchLength),fill=colors[2], bins=50)+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=branchLength),fill=colors[6], bins=50)#+
  #geom_histogram(data=subset(combined, outgroup==2),aes(x=branchLength),fill=colors[4], bins=50)
branchLengths

introProbs <- ggplot() +
  geom_histogram(data=subset(combined, outgroup==3),aes(x=introProb),fill=colors[2])+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=introProb),fill=colors[6])+
  geom_histogram(data=subset(combined, outgroup==2),aes(x=introProb),fill=colors[4])
introProbs

BLvsIP <- ggplot(data=combined) +
  geom_point(aes(x=branchLength,y=introProb, col=as.factor(outgroup)))+
  scale_color_manual(values=c(colors[6],colors[4],colors[2]))
BLvsIP

BLvsFd <- ggplot(data=combined) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.7)+
  scale_color_manual(values=c(colors[6],colors[4],colors[2]))
BLvsFd

BLvsFd_intro <- ggplot(data=subset(combined,outgroup==1)) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.7)+
  scale_color_manual(values=c(colors[6]))
BLvsFd_intro
ggsave(file="HeliconiusComparativeGenomics/FdComparison/5KB_hq_branchLengthVsFd_intro.png",plot=BLvsFd_intro, height=10,width=10)


BLvsFd_species <- ggplot(data=subset(combined,outgroup==3)) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.7)+
  scale_color_manual(values=c(colors[2]))
BLvsFd_species
ggsave(file="HeliconiusComparativeGenomics/FdComparison/5KB_hq_branchLengthVsFd_species.png",plot=BLvsFd_species, height=10,width=10)

BLvsFd_ILS <- ggplot(data=subset(combined,outgroup==2)) +
  geom_point(aes(x=branchLength,y=fd, col=as.factor(outgroup)), alpha=0.7)+
  scale_color_manual(values=c(colors[4]))
BLvsFd_ILS
ggsave(file="HeliconiusComparativeGenomics/FdComparison/5KB_hq_branchLengthVsFd_ILS.png",plot=BLvsFd_ILS, height=10,width=10)


IPvsFd <- ggplot(data=combined) +
  geom_point(aes(x=introProb,y=fd, col=as.factor(outgroup)))+
  scale_color_manual(values=c(colors[6],colors[4],colors[2]))
IPvsFd
ggsave(file="HeliconiusComparativeGenomics/FdComparison/5KB_hq_introProbVsFd_ILS.png",plot=IPvsFd, height=10,width=10)


##### Plot mean Fd by chromosome size ########

meanFd <- c()
for (chr in HeratoGenome$chromosome){
  meanFd <- c(meanFd,mean(abs(treeRate[which(treeRate$chrom.x == chr),]$fd)))
}
HeratoGenome$meanFd <- meanFd
HeratoGenome$lengthMB <- HeratoGenome$end/1000000
length_rho <- cor.test(~meanFd+lengthMB, data=subset(HeratoGenome, chr !=21), method="spearman")
length_lm <- lm(meanFd~lengthMB, data=subset(HeratoGenome, chr !=21))
colors=brewer.pal(8,"Paired")
plot <- ggplot(data=HeratoGenome) +
  geom_point(aes(y=meanFd,x=lengthMB),col=colors[2]) +
  geom_text_repel(aes(label=chromosome,y=rep(max(meanFd+.005),21),x=lengthMB),direction="x",nudge_y=.01)+
  geom_abline(slope=length_lm$coefficients[2],intercept=length_lm$coefficients[1],col=colors[2])+
  # #ylim(0,0.70)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Length (MB)",y=expression(italic("f"['d'])))
plot




##### plot Fd by local recombination rate #######

#unbinned
fdModel <- lm(fd~LD, data=treeRate)
lmIntercept=fdModel$coefficients[1]
lmSlope=fdModel$coefficients[2]
confFrame <- data.frame(LD=treeRate$LD)
confInterval <- as.data.frame(predict(fdModel, confFrame, interval="confidence"))
fdVsRate <- ggplot(data=treeRate)+
  geom_point(aes(x=LD,y=fd), alpha=0.2)+
  geom_line(aes(x=confFrame$LD, y=confInterval$lwr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$upr), lty=3, col="orange")+
  geom_line(aes(x=confFrame$LD, y=confInterval$fit), lty=1, col="orange")+
  geom_smooth(aes(x=LD,y=fd))+
  labs(x="recombination rate",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
fdVsRate


##### separate into bins ########
rateKW <- kruskal(treeRate$fd, treeRate$rateQuint, group=TRUE, p.adj="bonferroni")
ratedf <- data.frame(tree=rownames(rateKW$groups), rateQuint=rateKW$groups$groups)

FdByRate <-   ggplot(data=treeRate)+
  geom_boxplot(aes(y=fd, x=rateQuint),outlier.shape = NA)+
  geom_jitter(aes(y=fd,x=rateQuint), alpha=.05)+
  geom_text(data=ratedf, aes(x=tree, y=0.8, label=rateQuint), size=8)+
  labs(x="recombination rate quintile",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
FdByRate

######### Fd by tree ########
treeSubsetFrame <- data.frame()
for (treeType in unique(bestTrees$tree)){
  treeSub <- bestTrees[sample(which(bestTrees$tree==treeType),300, replace = F),]
  treeSubsetFrame <- rbind(treeSubsetFrame,treeSub )
}

treeKW <- kruskal(treeSubsetFrame$fd, treeSubsetFrame$tree, group=TRUE, p.adj="bonferroni")
treedf <- data.frame(tree=rownames(treeKW$groups), treeGroup=treeKW$groups$groups)
treedf$tree <- factor(treedf$tree, levels=fiveEveryFive)

FdByTree <-   ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=fd, x=tree),outlier.shape = NA)+
  geom_jitter(aes(y=fd,x=tree), alpha=.05)+
  scale_x_discrete(labels=c("Tree 1","Tree 2","Tree 3","Tree 4","Tree 5","Tree 6","Tree 7","Tree 8"))+
  geom_text(data=treedf, aes(x=tree, y=0.6, label=treeGroup), size=8)+
  labs(x="",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
FdByTree


####### D by tree #########
DByTree <-   ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=D, x=tree),outlier.shape = NA)+
  geom_jitter(aes(y=D,x=tree), alpha=.05)+
  scale_x_discrete(labels=c("Tree 1","Tree 2","Tree 3","Tree 4","Tree 5","Tree 6","Tree 7","Tree 8"))+
  #geom_text(data=treedf, aes(x=tree, y=0.6, label=treeGroup), size=8)+
  labs(x="",y=expression(italic("D")))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12))
DByTree

########## fd by miyagi call #####

callSubsetFrame <- data.frame()
for (callType in unique(treeRate$isIntro)){
  callSub <- treeRate[sample(which(treeRate$isIntro==callType),400, replace = F),]
  callSubsetFrame <- rbind(callSubsetFrame,callSub )
}

callKW <- kruskal(callSubsetFrame$fd, callSubsetFrame$isIntro, group=TRUE, p.adj="bonferroni")
calldf <- data.frame(call=rownames(callKW$groups), callGroup=callKW$groups$groups)

FdByCall <-   ggplot(data=callSubsetFrame)+
  geom_boxplot(aes(y=fd, x=isIntro),outlier.shape = NA)+
  geom_jitter(aes(y=fd,x=isIntro), alpha=.05)+
  geom_text(data=calldf, aes(x=call, y=0.8, label=callGroup), size=8)+
  labs(x="",y=expression(italic("f"['d'])))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=12),axis.text.x = element_text(angle = 45, hjust = 1))
FdByCall

####### make figure #######

ng <- nullGrob()
g <- arrangeGrob(FdByRate, FdByTree, 
                 FdByCall, ng)
ggsave(filename = "fd_recomb_miyagi5KBevery50_kwtest.pdf",g, width = 10, height=10)
