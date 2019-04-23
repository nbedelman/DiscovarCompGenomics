#Script to generate Supplementary Figure XX

############## Import Libraries ###############
library(devtools)
library("karyoploteR")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plotly)
library(tidyr)
library(knitr)
library(gridExtra)
library(grid)
library(agricolae)
library(ape)
library(phangorn)
library(foreach)
library(doParallel)



setwd("/Users/nbedelman/Dropbox/HeliconiusComparativeGenomics/")

########### Define Data Files #############
#10KB abutting windows - Figure SXX
STATS="datasets/eratoClade_10KBAbutting_HeraRef_windowStats.txt"
STATS2="Figure3/eratoClade_HeraRef_fullAlignedSites.windowStats.txt"
TREES="Figure3/eratoClade_HeraRef_allSites.trees.tsv"
TREES2="Figure3/eratoClade_HeraRef_fullAlignedSites.trees.tsv"
out="eratoClade_10KBAbutting_HeraRef"

############# Define helper functions ############
getChrom <- function (x) {
  chrom <- strsplit(x,"_")[[1]][1]
  return(chrom)
}

getEnd <- function (x) {
  end <- strsplit(x,"_")[[1]][2]
  return(end)
}


############## Read in data #########
#trees
trees <- read.csv(TREES, header=FALSE, col.names = c("block_id","tree"), sep=" ")
trees$block_id <- as.character(trees$block_id)

trees2 <- read.delim(TREES2, header=FALSE, col.names = c("block_id","tree"), sep=" ")
trees2$block_id <- as.character(trees2$block_id)


#stats
stats <- read.delim(STATS)
stats2=read.delim(STATS2)
stats=inner_join(stats,stats2, by="block_id")
stats$block_id <- as.character(stats$block_id)
stats$mean_bootstrap_support.x=as.numeric(as.character(stats$mean_bootstrap_support.x))
stats$mean_bootstrap_support.y=as.numeric(as.character(stats$mean_bootstrap_support.y))


############### Tidy Data and make overview graphs #########
allData <- inner_join(trees,stats) 
allData <- inner_join(allData, trees2, by="block_id")
#allData <- stats
allData <- allData %>% drop_na
allData$phi_p.x <- as.numeric(as.character(allData$phi_p.x))
allData$transformedP.x <- -log10(allData$phi_p.x)
allData$phi_p.y <- as.numeric(as.character(allData$phi_p.y))
allData$transformedP.y <- -log10(allData$phi_p.y)
allData$tree.x <- as.character(allData$tree.x)
allData$tree.y <- as.character(allData$tree.y)


##Compare Trees #####
distances <- foreach(rw=seq(1,nrow(allData)),.combine=c)%do%{
  tree1=read.tree(text=allData$tree.x[rw])
  tree2=read.tree(text=allData$tree.y[rw])
  dist=RF.dist(as.phylo(tree1), as.phylo(tree2))
}


allData$RF <- distances

############ Compare datasets ################
bootCompare <- ggplot(data=allData)+
geom_point(aes(x=allData$mean_bootstrap_support.x, y=allData$mean_bootstrap_support.y), alpha=.2)+
  geom_abline(aes(slope=1, intercept=0))
bootCompare

phiCompare <- ggplot(data=allData)+
  geom_point(aes(x=allData$transformedP.x, y=allData$transformedP.y), alpha=.2)+
  geom_abline(aes(slope=1, intercept=0))
phiCompare

RFvsBoot <- ggplot(data=allData)+
  geom_point(aes(x=allData$mean_bootstrap_support.x, y=allData$RF), alpha=.2)
RFvsBoot

allData <- subset(allData, mean_bootstrap_support.x >=80 & fullLength.x>=2000)
hist(allData$RF)
allData <- subset(allData, fullLength>100)

lengths <- ggplot(data=allData)+
  geom_histogram(aes(x=fullLength.x))+
  labs(x="Alignment length", y="Number of windows")

pctMissing <- ggplot(data=allData)+
  geom_histogram(aes(x=p_missing.x))+
  labs(x="Fraction missing", y="Number of windows")

piSites <- ggplot(data=allData)+
  geom_histogram(aes(x=n_pi_sites.x))+
  labs(x="Phylogenetically informative sites", y="Number of windows")

bootstrap <- ggplot(data=allData)+
  geom_histogram(aes(x=mean_bootstrap_support))+
  labs(x="Mean bootstrap support", y="Number of windows")

breakpoint <- ggplot(data=allData)+
  geom_histogram(aes(x=transformedP))+
  labs(x="Probability of recombination within block", y="Number of windows")

completeSites <- ggplot(data=allData)+
  geom_histogram(aes(x=fullAlignedSites))+
  labs(x="Fully aligned sites", y="Number of windows")

ng <- nullGrob()

grid.arrange(lengths, pctMissing, 
             piSites,bootstrap,
             breakpoint,completeSites)

# grid.arrange(lengths, pctMissing, 
#              piSites,
#              breakpoint,ng)

g <- arrangeGrob(lengths, pctMissing, 
                 piSites,bootstrap,
                 breakpoint,completeSites)
# 
# g <- arrangeGrob(lengths, pctMissing, 
#                  piSites,
#                  breakpoint,ng)
ggsave(filename = paste0(out,"individualStats.pdf"),g, width = 10, height=10)
#Pairwise Graphs

oneOne <- ggplot(data=allData)+
  geom_point(aes(x=fullLength,y=n_pi_sites),alpha=.05)+
  labs(y="Phylogenetically informative sites",x="")

oneTwo<- ggplot(data=allData)+
  geom_point(aes(x=fullLength,y=p_missing),alpha=.05)+
  labs(y="Fraction missing",x="")

oneThree<- ggplot(data=allData)+
  geom_point(aes(x=fullLength,y=transformedP),alpha=.05)+
  labs(y="Probability of recombination within block",x="")

oneFour <- ggplot(data=allData)+
  geom_point(aes(x=fullLength,y=mean_bootstrap_support),alpha=.05)+
  labs(y="Mean bootstrap support",x="Alignment length")

twoTwo <- ggplot(data=allData)+
  geom_point(aes(x=n_pi_sites,y=p_missing),alpha=.05)+
  labs(y="",x="")

twoThree <- ggplot(data=allData)+
  geom_point(aes(x=n_pi_sites,y=transformedP),alpha=.05)+
  labs(y="",x="")

twoFour <- ggplot(data=allData)+
  geom_point(aes(x=n_pi_sites,y=mean_bootstrap_support),alpha=.05)+
  labs(y="",x="Phylogenetically informative sites")

threeThree <- ggplot(data=allData)+
  geom_point(aes(x=p_missing,y=transformedP),alpha=.05)+
  labs(y="",x="")

threeFour <- ggplot(data=allData)+
  geom_point(aes(x=p_missing,y=mean_bootstrap_support),alpha=.05)+
  labs(y="",x="Fraction missing")

fourFour <- ggplot(data=allData)+
  geom_point(aes(x=transformedP,y=mean_bootstrap_support),alpha=.05)+
  labs(y="",x="Probability of recombination within block")

ng <- nullGrob()
grid.arrange(oneOne, ng, ng,ng,        
             oneTwo, twoTwo,ng,ng,
             oneThree,twoThree,threeThree,ng,
             oneFour,twoFour,threeFour,fourFour)

g <- arrangeGrob(oneOne, ng, ng,ng,        
                 oneTwo, twoTwo,ng,ng,
                 oneThree,twoThree,threeThree,ng,
                 oneFour,twoFour,threeFour,fourFour)
ggsave(filename = paste0(out,"allByAll.jpg"),g, width = 13, height=13)


######## analyze comparisons between trees ###########
colors=brewer.pal(8,"Paired")
treeLabels <- c("tree0","tree1","tree7","tree2","tree5","tree3","tree4","tree24")
colorOrder=c(colors[6],colors[2],colors[3],colors[4],colors[8],colors[1],colors[7],colors[5])
treeLabels=c("tree1", "tree2", "tree3", "tree4", "tree5", "tree6", "tree7", "tree8")
topTrees <- subset(allData, tree %in% treeLabels& fullLength>2000)# & mean_bootstrap_support>=75 & p_missing <.4)
topTrees$tree <- factor(topTrees$tree, levels <- treeLabels)

treeSubsetFrame <- data.frame()
for (treeType in treeLabels){
  treeSub <- topTrees[sample(which(topTrees$tree==treeType),500, replace = F),]
  treeSubsetFrame <- rbind(treeSubsetFrame,treeSub )
}

lengthKW <- kruskal(treeSubsetFrame$fullLength, treeSubsetFrame$tree, group=TRUE, p.adj="bonferroni")
lengthdf <- data.frame(tree=rownames(lengthKW$groups), lengthGroup=lengthKW$groups$groups)

piSitesKW <- kruskal(treeSubsetFrame$n_pi_sites, treeSubsetFrame$tree, group=TRUE, p.adj="bonferroni")
piSitesdf <- data.frame(tree=rownames(piSitesKW$groups), piSitesGroup=piSitesKW$groups$groups)

missingKW <- kruskal(treeSubsetFrame$p_missing, treeSubsetFrame$tree, group=TRUE, p.adj="bonferroni")
missingdf <- data.frame(tree=rownames(missingKW$groups), missingGroup=missingKW$groups$groups)

bootKW <- kruskal(treeSubsetFrame$mean_bootstrap_support, treeSubsetFrame$tree, group=TRUE, p.adj="bonferroni")
bootdf <- data.frame(tree=rownames(bootKW$groups), bootGroup=bootKW$groups$groups)

recombKW <- kruskal(treeSubsetFrame$transformedP, treeSubsetFrame$tree, group=TRUE, p.adj="bonferroni")
recombdf <- data.frame(tree=rownames(recombKW$groups), recombGroup=recombKW$groups$groups)

KWdf <- inner_join(lengthdf,piSitesdf, by="tree") %>% inner_join(., missingdf, by="tree") %>% inner_join(., bootdf, by="tree") %>% inner_join(., recombdf, by="tree")
KWdf$tree <- factor(KWdf$tree,levels= treeLabels)

bOneOne <-  ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=fullLength,x=tree, color=tree), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=fullLength,x=tree), alpha=.05)+
  scale_color_manual(values =colorOrder)+
  scale_x_discrete(labels=treeLabels)+
  geom_text(data=KWdf, aes(x=tree, y=10500, label=lengthGroup), size=8)+
  labs(y="Length",x="")+
  guides(col=FALSE)

bOneTwo <-  ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=n_pi_sites,x=tree, color=tree), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=n_pi_sites,x=tree), alpha=.05)+
  scale_color_manual(values =colorOrder)+
  scale_x_discrete(labels=treeLabels)+
  geom_text(data=KWdf, aes(x=tree, y=1.05*max(topTrees$n_pi_sites), label=piSitesGroup), size=8)+
  labs(y="Phylogenetic informative sites",x="")+
  guides(col=FALSE)

bOneThree <-  ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=transformedP,x=tree, color=tree), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=transformedP,x=tree), alpha=.05)+
  scale_color_manual(values =colorOrder)+
  scale_x_discrete(labels=treeLabels)+
  geom_text(data=KWdf, aes(x=tree, y=30, label=recombGroup), size=8)+
  labs(y="Probability of recombination",x="")+
  guides(col=FALSE)

bTwoOne <-  ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=mean_bootstrap_support,x=tree, color=tree), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=mean_bootstrap_support,x=tree), alpha=.05)+
  scale_color_manual(values =colorOrder)+
  scale_x_discrete(labels=treeLabels)+
  geom_text(data=KWdf, aes(x=tree, y=1.05*max(treeSubsetFrame$mean_bootstrap_support), label=bootGroup), size=8)+
  labs(y="Bootstrap value",x="")+
  guides(col=FALSE)


bTwoTwo <-  ggplot(data=treeSubsetFrame)+
  geom_boxplot(aes(y=p_missing,x=tree, color=tree), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=p_missing,x=tree), alpha=.05)+
  scale_color_manual(values =colorOrder)+
  scale_x_discrete(labels=treeLabels)+
  geom_text(data=KWdf, aes(x=tree, y=1.05*max(treeSubsetFrame$p_missing), label=missingGroup), size=8)+
  labs(y="Percent missing data",x="")+
  guides(col=FALSE)



grid.arrange(bOneOne, bOneTwo, 
             bTwoOne,bTwoTwo,
             bOneThree,ng)
g <- arrangeGrob(bOneOne, bOneTwo, 
                 bTwoOne,bTwoTwo,
                 bOneThree,ng)
ggsave(filename = paste0(out,"statsByTree.jpg"),g, width = 10, height=10)

## Analyze comparisons between recombination rates ## 
######## analyze comparisons between trees ###########
lengthKW <- kruskal(allData$fullLength, allData$rateQuint, group=TRUE, p.adj="bonferroni")
lengthdf <- data.frame(rateQuint=rownames(lengthKW$groups), lengthGroup=lengthKW$groups$groups)

piSitesKW <- kruskal(allData$n_pi_sites, allData$rateQuint, group=TRUE, p.adj="bonferroni")
piSitesdf <- data.frame(rateQuint=rownames(piSitesKW$groups), piSitesGroup=piSitesKW$groups$groups)

missingKW <- kruskal(allData$p_missing, allData$rateQuint, group=TRUE, p.adj="bonferroni")
missingdf <- data.frame(rateQuint=rownames(missingKW$groups), missingGroup=missingKW$groups$groups)

bootKW <- kruskal(allData$mean_bootstrap_support, allData$rateQuint, group=TRUE, p.adj="bonferroni")
bootdf <- data.frame(rateQuint=rownames(bootKW$groups), bootGroup=bootKW$groups$groups)

recombKW <- kruskal(allData$transformedP, allData$rateQuint, group=TRUE, p.adj="bonferroni")
recombdf <- data.frame(rateQuint=rownames(recombKW$groups), recombGroup=recombKW$groups$groups)

KWdf <- inner_join(lengthdf,piSitesdf, by="rateQuint") %>% inner_join(., missingdf, by="rateQuint") %>% inner_join(., bootdf, by="rateQuint") %>% inner_join(., recombdf, by="rateQuint")
KWdf$rateQuint <- factor(KWdf$rateQuint,levels= c("1","2","3","4","5"))

rOneOne <-  ggplot(data=allData)+
  geom_boxplot(aes(y=fullLength,x=rateQuint), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=fullLength,x=rateQuint),col="grey", alpha=.05)+
  geom_text(data=KWdf, aes(x=rateQuint, y=10500, label=lengthGroup), size=8)+
  labs(y="Length",x="")+
  guides(col=FALSE)

rOneTwo <-  ggplot(data=allData)+
  geom_boxplot(aes(y=n_pi_sites,x=rateQuint), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=n_pi_sites,x=rateQuint), col="grey",alpha=.05)+
  geom_text(data=KWdf, aes(x=rateQuint, y=1.05*max(allData$n_pi_sites), label=piSitesGroup), size=8)+
  labs(y="Phylogenetic informative sites",x="")+
  guides(col=FALSE)

rOneThree <-  ggplot(data=allData)+
  geom_boxplot(aes(y=transformedP,x=rateQuint), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=transformedP,x=rateQuint), col="grey",alpha=.05)+
  geom_text(data=KWdf, aes(x=rateQuint, y=30, label=recombGroup), size=8)+
  labs(y="Probability of recombination",x="")+
  guides(col=FALSE)

rTwoOne <-  ggplot(data=allData)+
  geom_boxplot(aes(y=mean_bootstrap_support,x=rateQuint), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=mean_bootstrap_support,x=rateQuint), col="grey",alpha=.05)+
  geom_text(data=KWdf, aes(x=rateQuint, y=1.05*max(allData$mean_bootstrap_support), label=bootGroup), size=8)+
  labs(y="Bootstrap value",x="")+
  guides(col=FALSE)


rTwoTwo <-  ggplot(data=allData)+
  geom_boxplot(aes(y=p_missing,x=rateQuint), lwd=2,outlier.shape = NA)+
  geom_jitter(aes(y=p_missing,x=rateQuint), col="grey",alpha=.05)+
  geom_text(data=KWdf, aes(x=rateQuint, y=1.05*max(allData$p_missing), label=missingGroup), size=8)+
  labs(y="Percent missing data",x="")+
  guides(col=FALSE)



grid.arrange(rOneOne, rOneTwo, 
             rTwoOne,rTwoTwo,
             rOneThree,ng)
g <- arrangeGrob(rOneOne, rOneTwo, 
                 rTwoOne,rTwoTwo,
                 rOneThree,ng)
ggsave(filename = paste0(out,"statsByRecombRate.jpg"),g, width = 10, height=10)

####### Plot data to H. erato ###########
ABGranges <- makeGRangesFromDataFrame(allData,keep.extra.columns = TRUE)
ABGranges <- sortSeqlevels(ABGranges)
ABGranges <- sort(ABGranges)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight=0

kp <- plotKaryotype(genome=HeratoGenomeGrange,plot.type=1,cex=2,plot.params = pp)


#Trees are assigned arbitrary numbers in the filtering process, so need to assign the proper colors separately for each dataset
############# Plot trees for 10KB blocks ############
#window trees
kpDataBackground(kp, data.panel = 1,r0=0, r1=.75, col="black")
kpRect(kp, subset(ABGranges,tree=="tree1"),y0=0,y1=1,r0=0, r1=.75, col=colors[2], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree0"),y0=0,y1=1,r0=0,r1=.75, col=colors[6], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree7"),y0=0,y1=1,r0=0, r1=.75, col=colors[3], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree2"),y0=0,y1=1,r0=0, r1=.75, col=colors[4], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree4"),y0=0,y1=1,r0=0, r1=.75, col=colors[7], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree5"),y0=0,y1=1,r0=0, r1=.75, col=colors[8], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree3"),y0=0,y1=1,r0=0, r1=.75, col=colors[1], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree24"),y0=0,y1=1,r0=0, r1=.75, col=colors[5], border=NA)
kpRect(kp, subset(ABGranges,! tree %in% c("tree1","tree0","tree7","tree2","tree4","tree5","tree3","tree24")),
       y0=0,y1=1,r0=0, r1=0.75, col="grey", border=NA)


