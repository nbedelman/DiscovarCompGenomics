#### compare miyagi output to Fd in simulations

library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(agricolae)

setwd("~/Dropbox")

FD="HeliconiusComparativeGenomics/FdComparison/fd_sim1Output.txt"
MIYAGI <- "HeliconiusComparativeGenomics/FdComparison/eratoClade_HeraRef_5KB_pRecombLe.05.miyagi.csv"

########## Define Functions #############

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

getSegment <- function (x) {
  segment <- strsplit(x,"_")[[1]][2]
  return(segment)
}


colors=brewer.pal(8,"Paired")


############ Read and format data #############

fd <- read.csv(FD)
fd$segment <- as.character(unlist(lapply(as.character(fd$scaffold),getSegment)))
fd$segment <- as.numeric(fd$segment)-1

miyagi <- read.csv(MIYAGI, header=F, col.names = c("segment","outgroup","branchLength","introProb","blank"))
miyagi <- select(miyagi, -"blank")


combined <- inner_join(miyagi,fd)
combined <- miyagi
############ plot basic graphs ##########
branchLengths <- ggplot()+
  geom_histogram(data=subset(combined, outgroup==3),aes(x=branchLength), fill=colors[2],bins =20)+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=branchLength), fill=colors[6], bins=20)+
  geom_histogram(data=subset(combined, outgroup==2),aes(x=branchLength), fill=colors[4], bins=20)+
  scale_x_continuous(limits=c(0,0.01))
branchLengths

introProbs <- ggplot()+
  geom_histogram(data=subset(combined, outgroup==3),aes(x=introProb), fill=colors[2], alpha=0.5)+
  geom_histogram(data=subset(combined, outgroup==1),aes(x=introProb), fill=colors[6], alpha=0.75)+
  geom_histogram(data=subset(combined, outgroup==2),aes(x=introProb), fill=colors[4])
introProbs

introVSBranch <- ggplot()+
  geom_point(data=subset(combined, outgroup==3),aes(x=branchLength,y=introProb), col=colors[2])+
  geom_point(data=subset(combined, outgroup==1),aes(x=branchLength,y=introProb), col=colors[6])+
  geom_point(data=subset(combined, outgroup==2),aes(x=branchLength,y=introProb), col=colors[4])+
  scale_x_continuous(limits=c(0,0.1))
introVSBranch


introVSFd <- ggplot()+
  geom_boxplot(data=combined,aes(x=outgroup, group=outgroup, y=fd, col=as.factor(outgroup))) +
  scale_color_manual(values=c(colors[4],colors[6],colors[2]))
introVSFd

introVSD <- ggplot()+
  geom_boxplot(data=combined,aes(x=outgroup, group=outgroup, y=D, col=as.factor(outgroup))) +
  scale_color_manual(values=c(colors[4],colors[6],colors[2]))
introVSD


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
