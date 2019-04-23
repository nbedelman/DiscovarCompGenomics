library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(dplyr)
library(devtools)


setwd("/Users/nbedelman/dropbox/HeliconiusComparativeGenomics/")

TREES="datasets/eratoClade.10KBAbutting.iqTree.sorted.csv"
STATS="phylogenies_stats/10KB_eratoClade/eratoClade_10KBAbutting_HeraRef_windowStats.txt"
RecRate="datasets/HeratoWindows.10KBAbutting.recRates.CDSoverlaps.bed"
GENOMEFILE="datasets/Herato.transitions_fullChroms.bed"
out="Figure3/eratoClade_10KBAbutting_ge2kb_le.2miss_"

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

addTreeToFrame <- function(treeName,frame){
  chroms <- HeratoGenome$chromosome
  treeNums <- c()
  for (chr in chroms){
  treeNums <- c(treeNums,length(which(trees$tree==treeName & trees$chrom ==chr))/length(which(trees$chrom ==chr)))
  }
frame[[treeName]] <- treeNums
return(frame)
}

getChromLength <- function(chr){
  return(HeratoGenome[which(HeratoGenome$chromosome==chr),]$end)
}

addToCountsFrame <- function(t,frame,numBreaks){
  frame[[t]] <- hist(subset(trees,tree==t)$distFromCenter, breaks = numBreaks,plot=F)$counts
  return(frame)
}

addToFrame <- function(t,frame,numBreaks,variable){
  frame[[t]] <- hist(subset(trees,tree==t)[[variable]], breaks = numBreaks,plot=F)$counts
  return(frame)
}
########## read and format data #############

HeratoGenome <- read.delim(GENOMEFILE, col.names=c("seqnames","start","end","chromosome"), header=F)
HeratoGenome$chr <-  unlist(strsplit(as.character(HeratoGenome$chromosome),"chr"))[c(FALSE,TRUE)]

LDRates <- "../heliconius/MANUSCRIPT/Sections/recombinationRate/LDhelmet_erato.txt"
LDMap <- read.delim(LDRates)
LDMap$chrom <- unlist(strsplit(as.character(LDMap$scaffold),"Herato"))[c(FALSE,TRUE)]
LDMap$chrom <- substr(LDMap$chrom,rep(1,length(LDMap$chrom)),rep(2,length(LDMap$chrom)))
LDRate <- c()
for (i in unique(LDMap$chrom)){
  chromSub <- subset(LDMap,chrom==i)
  LDRate <- c(LDRate,mean(chromSub$meanS/1000000, na.rm=T))
}

mapFrame <- data.frame(chr=as.character(seq(1,21)), mapLength_cross=c(49,38,47,48,44,44,42,40,42,47,45,50,44,39,52,47,36,52,47,52,41), recRate_LD=LDRate)


HeratoGenome <- inner_join(HeratoGenome,mapFrame)
HeratoGenome$recRate_cross=HeratoGenome$mapLength_cross/HeratoGenome$end
HeratoGenome$recRate_cross_MB <- HeratoGenome$recRate_cross*1000000
HeratoGenome$mapLength_LD <- HeratoGenome$recRate_LD*HeratoGenome$end
HeratoGenome$recRate_LD_MB <- HeratoGenome$recRate_LD*1000000



trees <- read.csv(TREES, header=FALSE, col.names = c("segment","tree"))
trees$segment <- as.character(trees$segment)
trees$chrom <- as.character(unlist(lapply(trees$segment,getChrom)))
trees$end <- as.numeric(unlist(lapply(trees$segment,getEnd)))*10000
trees$start <- trees$end-9999


stats <- read.delim(STATS)
stats$segment <- as.character(stats$block_id)
stats$mean_bootstrap_support=as.numeric(as.character(stats$mean_bootstrap_support))

colors=brewer.pal(8,"Paired")


########## Filter Data ###########
trees <- inner_join(trees,stats)
trees <- subset(trees, fullLength>=2000 & p_missing <.2)

########## tree fraction to chromosome size; figure 3A, S7.3  #############
#non full-aligned: 0,1,7,2,4,24,5,3
#full-aligned: 2,3,6,0,1,15,5,10
topTrees <- c("tree0","tree1","tree7","tree2","tree4","tree24","tree5","tree3")
topTwo <- c("tree0","tree1")
colorOrder <- c(colors[6],colors[2],colors[3],colors[4],colors[7],colors[5],colors[8],colors[1])

lengthFrame <- data.frame(chroms=HeratoGenome$chromosome,length=HeratoGenome$end/1000000, crossLength=HeratoGenome$recRate_cross_MB)
lengthFrame <- addTreeToFrame("tree0",lengthFrame) # tel-hsa sister (red)
lengthFrame <- addTreeToFrame("tree1",lengthFrame) # species tree (blue)
lengthFrame <- addTreeToFrame("tree7",lengthFrame) # light green
lengthFrame <- addTreeToFrame("tree2",lengthFrame) # dark green
lengthFrame <- addTreeToFrame("tree4",lengthFrame) # light orange
lengthFrame <- addTreeToFrame("tree24",lengthFrame) # light red
lengthFrame <- addTreeToFrame("tree5",lengthFrame) # dark orange
lengthFrame <- addTreeToFrame("tree3",lengthFrame) # light blue


lengthFrameSub <- subset(lengthFrame, ! chroms %in% c("chr21","chr0"))
########## find correlations #############
spearmanRho <- data.frame("tree"=vector(),"rho"=vector(),"p"=vector(),"slope"=vector(),"intercept"=vector())
for(column in topTrees){
  fmla <- as.formula(paste0(column,"~ length"))
  length_fit <- cor.test(fmla, data=lengthFrameSub, method="spearman")
  corrs <- rbind(corrs,data.frame("tree"=column,"r2"=round(summary(length_fit)$adj.r.squared,3),"p"=round(lmp(length_fit),4),"slope"=length_fit$coefficients[2],"intercept"=length_fit$coefficients[1]))
}


corrs <- data.frame("tree"=vector(),"r2"=vector(),"p"=vector(),"slope"=vector(),"intercept"=vector())
for(column in topTrees){
  fmla <- as.formula(paste0(column,"~ length"))
  length_fit <- lm(fmla, data=lengthFrameSub)
  corrs <- rbind(corrs,data.frame("tree"=column,"r2"=round(summary(length_fit)$adj.r.squared,3),"p"=round(lmp(length_fit),4),"slope"=length_fit$coefficients[2],"intercept"=length_fit$coefficients[1]))
}


########## plot chromosome proportions #############
###Top 2 (for figure 3)
colors=brewer.pal(8,"Paired")
plot <- ggplot(data=lengthFrame) +
  geom_point(aes(y=tree0,x=length),col=colors[6]) +
  geom_point(aes(y=tree1,x=length),col=colors[2])+
  geom_abline(slope=corrs$slope[1],intercept=corrs$intercept[1],col=colors[6])+
  geom_abline(slope=corrs$slope[2],intercept=corrs$intercept[2],col=colors[2])+
  geom_text_repel(aes(label=chroms,y=0.65,x=length),direction="x",nudge_y=.03)+
  ylim(0,0.70)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Length (MB)",y="Fraction of Chromosome")
plot
ggsave(plot, filename = paste0(out,"treesByChromSize_topTwo.pdf"), width = 10, height=10)


##all 

plot <- ggplot(data=lengthFrame) +
  geom_point(aes(y=tree0,x=length),col=colors[6]) +
  geom_point(aes(y=tree1,x=length),col=colors[2]) +
  geom_point(aes(y=tree7,x=length),col=colors[3]) +
  geom_point(aes(y=tree2,x=length),col=colors[4]) +
  geom_point(aes(y=tree4,x=length),col=colors[7]) +
  geom_point(aes(y=tree24,x=length),col=colors[5]) +
  geom_point(aes(y=tree5,x=length),col=colors[8]) +
  geom_point(aes(y=tree3,x=length),col=colors[1]) +
  geom_abline(slope=corrs$slope[1],intercept=corrs$intercept[1],col=colors[6])+
  geom_abline(slope=corrs$slope[2],intercept=corrs$intercept[2],col=colors[2])+
  geom_abline(slope=corrs$slope[3],intercept=corrs$intercept[3],col=colors[3])+
  geom_abline(slope=corrs$slope[4],intercept=corrs$intercept[4],col=colors[4])+
  geom_abline(slope=corrs$slope[5],intercept=corrs$intercept[5],col=colors[7])+
  geom_abline(slope=corrs$slope[6],intercept=corrs$intercept[6],col=colors[5])+
  geom_abline(slope=corrs$slope[7],intercept=corrs$intercept[7],col=colors[8])+
  geom_abline(slope=corrs$slope[8],intercept=corrs$intercept[8],col=colors[1])+
  geom_text_repel(aes(label=chroms,y=0.65,x=length),direction="x",nudge_y=.03)+
  ylim(0,0.70)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Length (MB)",y="Fraction of Chromosome")
plot

ggsave(plot, filename = paste0(out,"treesByChromSize_all.pdf"), width = 10, height=10)

######### plot chromosomal positions; figure 3B, S7.4, S7.5 ############
trees <- subset(trees,chrom != c("chr21", "chr0"))
trees$chromLengths <- unlist(lapply(trees$chrom,getChromLength))
trees$relPos <- ((trees$end+trees$start)/2)/trees$chromLengths
trees$distFromCenter <- abs(.5-trees$relPos)
trees <- subset(trees, distFromCenter>=0 & distFromCenter<=0.5)

interestingTrees <- subset(trees,tree %in% topTrees & ! chrom %in% c("chr21"))
interestingTrees$tree <- factor(interestingTrees$tree,levels=topTrees)

#### plot by window ###

breaks <- seq(0,.5,length.out = 11)
countsFrame <- data.frame(tree=vector(),position=vector(),chrom=vector(),prop=vector())
x <- hist(trees$distFromCenter, breaks = breaks,plot=F)$mids
for (chr in unique(interestingTrees$chrom)){
  #make a data frame for each chromosome
  chromTrees <- subset(trees,chrom == chr)
  total <- hist(chromTrees$distFromCenter, breaks = breaks,plot=F)$counts 
  for (t in topTrees){
    treeFrame <- subset(chromTrees, tree==t)
    treeCounts <- hist(treeFrame$distFromCenter, breaks = breaks,plot=F)$counts 
    entry <- data.frame(tree=rep(t,length(treeCounts)),position=x,chrom=rep(chr,length(treeCounts)),prop=treeCounts/total)
    countsFrame <- rbind(countsFrame,entry)
  }
}

countsFrame$position <- round(.5-countsFrame$position,2)
countsFrame$position <- as.factor(countsFrame$position)

###top two
plot <- ggplot(subset(countsFrame, tree %in% topTwo))+
  geom_boxplot(aes(x=position,y=prop,col=tree), show.legend=F) +
  labs(x="Distance from End (percent of chromosome length)",y="Fraction of Trees, by window")+
  scale_x_discrete(labels= c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50"))+
  theme(text=element_text(size=20))+
  scale_color_manual(values = colorOrder)
plot
ggsave(plot, filename = paste0(out,"treesByChromPosition_byWindow_topTwo.pdf"), width = 10, height=10)


plot <- ggplot(subset(countsFrame, tree %in% topTrees))+
  geom_boxplot(aes(x=position,y=prop,col=tree), show.legend=F) +
  labs(x="Distance from End (percent of chromosome length)",y="Fraction of Trees, by window")+
  scale_x_discrete(labels= c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40","40-45","45-50"))+
  theme(text=element_text(size=20))+
  scale_color_manual(values = colorOrder)
plot

ggsave(plot, filename = paste0(out,"treesByChromPosition_byWindow_all.pdf"), width = 10, height=10)

#### plot by toplogy 
breaks <- seq(0,.5,length.out = 11)
countsFrame2 <- data.frame(tree=vector(),position=vector(),chrom=vector(),prop=vector())
x <- hist(trees$distFromCenter, breaks = breaks,plot=F)$mids
for (chr in unique(interestingTrees$chrom)){
  #make a data frame for each chromosome
  chromTrees <- subset(trees,chrom == chr)
  for (t in topTrees){
    treeFrame <- subset(chromTrees, tree==t)
    total <- nrow(treeFrame)
    treeCounts <- hist(treeFrame$distFromCenter, breaks = breaks,plot=F)$counts 
    entry <- data.frame(tree=rep(t,length(treeCounts)),position=x,chrom=rep(chr,length(treeCounts)),prop=treeCounts/total)
    countsFrame2 <- rbind(countsFrame2,entry)
  }
}

##top two
countsFrame2$position <- round(.5-countsFrame2$position,2)
countsFrame2$position <- as.factor(countsFrame2$position)
plot <- ggplot(subset(countsFrame2, tree %in% topTwo))+
  geom_boxplot(aes(x=position,y=prop,col=tree), show.legend=F) +
  #xlim(c(0,1))+
  #facet_wrap(~position)+
  labs(x="Distance from End (fraction of chromosome length)",y="Fraction of Trees, by topology")+
  scale_x_discrete(labels= c("0-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50"))+
  scale_color_manual(values = colorOrder)+
  theme(text=element_text(size=20))
plot

ggsave(plot, filename = paste0(out,"treesByChromPosition_byTopolgy_topTwo.pdf"), width = 10, height=10)

###all trees
plot <- ggplot(subset(countsFrame2, tree %in% topTrees))+
  geom_boxplot(aes(x=position,y=prop,col=tree), show.legend=F) +
  #xlim(c(0,1))+
  #facet_wrap(~position)+
  labs(x="Distance from End (fraction of chromosome length)",y="Fraction of Trees, by topology")+
  scale_x_discrete(labels= c("0-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50"))+
  scale_color_manual(values = colorOrder)+
  theme(text=element_text(size=20))
plot

ggsave(plot, filename = paste0(out,"treesByChromPosition_byTopolgy_all.pdf"), width = 10, height=10)
############### Plot vs recombination Rate, Figure 3C, S7.6 ########

rateFile=read.delim(RecRate, col.names=c("scaffold","start","end","segment","dot","strand","loess2","loess3","loess4","LD","CDS"))

treeRate=inner_join(trees,rateFile, by="segment")
treeRate=treeRate[,c("segment","chrom","tree","loess2","loess3","loess4","LD","distFromCenter","CDS","chromLengths")]
treeRate$LD <- as.numeric(as.character(treeRate$LD))
treeRate$distFromCenter <- as.numeric(as.character(treeRate$distFromCenter))
treeRate$CDS <- as.numeric(as.character(treeRate$CDS))
treeRate=treeRate[complete.cases(treeRate),]
treeRate=subset(treeRate,LD>=0)
treeRate$CDSperCM <- (treeRate$CDS/treeRate$LD)/100

######### topology vs recombination rate #######
treeRate <- subset(treeRate,! chrom %in% c("chr21", "chr0"))
#breaks = quantile(treeRate$LD, probs = seq(0,1,.2)) #break by quintile
breaks=c(0,1,2,3,4,5,6,10)
rateFrame <- data.frame(tree=vector(),rate=vector(),chrom=vector(),prop=vector())
x <- hist(treeRate$LD, breaks = breaks,plot=F)$mids
nNum <- hist(treeRate$LD, breaks = breaks,plot=F)$counts
for (chr in unique(treeRate$chrom)){
  #make a data frame for each chromosome
  chromTrees <- subset(treeRate,chrom == chr)
  total <- hist(chromTrees$LD, breaks = breaks,plot=F)$counts 
  for (t in topTrees){
    treeFrame <- subset(chromTrees, tree==t)
    treeCounts <- hist(treeFrame$LD, breaks = breaks,plot=F)$counts 
    entry <- data.frame(tree=rep(t,length(treeCounts)),rate=x,chrom=rep(chr,length(treeCounts)),prop=treeCounts/total)
    rateFrame <- rbind(rateFrame,entry)
  }
}

rateFrame$rate <- round(rateFrame$rate,2)
rateFrame$rate<- as.factor(rateFrame$rate)
plot <- ggplot(subset(rateFrame, tree %in% topTwo))+
  geom_boxplot(aes(x=rate,y=prop,col=tree), show.legend=F, outlier.shape=NA) +
  labs(x="Recombination Rate (cM/MB)",y="Fraction of Trees")+
  scale_x_discrete(labels= c("1","2","3","4","5"))+ #quintiles
  scale_x_discrete(labels= c("0-1","1-2","2-3","3-4","4-5","5-6",">6"))+
  annotate("text",label=nNum,x=seq(1,length(breaks)-1),y=.8)+
  theme(text=element_text(size=20))+
  scale_color_manual(values = colorOrder)
plot

ggsave(plot, filename = paste0(out,"treesByRecombRate_topTwo.pdf"), width = 10, height=10)

plot <- ggplot(subset(rateFrame, tree %in% topTrees))+
  geom_boxplot(aes(x=rate,y=prop,col=tree), show.legend=F, outlier.shape=NA) +
  labs(x="Recombination Rate (cM/MB)",y="Fraction of Trees")+
  scale_x_discrete(labels= c("1","2","3","4","5"))+ #quintiles
  scale_x_discrete(labels= c("0-1","1-2","2-3","3-4","4-5","5-6",">6"))+
  annotate("text",label=nNum,x=seq(1,length(breaks)-1),y=.8)+
  theme(text=element_text(size=20))+
  scale_color_manual(values = colorOrder)
plot
ggsave(plot, filename = paste0(out,"treesByRecombRate_all.pdf"), width = 10, height=10)


for(bin in unique(rateFrame$rate)){
  print(bin)
  print(t.test(paired=T, x=subset(rateFrame, rate==bin & tree==topTwo[1])$prop,y=subset(rateFrame, rate==bin & tree==topTwo[2])$prop))
}


######### topology vs gene density #######
treeRate <- subset(treeRate,! chrom %in% c("chr21"))
# treeRate$logCDSdens <- log(treeRate$CDSperCM+1)
# breaks = quantile(treeRate$CDSperCM, probs = c(0,seq(.5,1,.1)))
treeRate$logCDSdens <- log(treeRate$CDS+1)
breaks = c(0,1,4,5,6,7,8,10)
cdsFrame <- data.frame(tree=vector(),cds=vector(),chrom=vector(),prop=vector())
x <- hist(treeRate$logCDSdens, breaks = breaks,plot=F)$mids
for (chr in unique(treeRate$chrom)){
  #make a data frame for each chromosome
  chromTrees <- subset(treeRate,chrom == chr)
  total <- hist(chromTrees$logCDSdens, breaks = breaks,plot=F)$counts 
  for (t in topTrees){
    treeFrame <- subset(chromTrees, tree==t)
    treeCounts <- hist(treeFrame$logCDSdens, breaks = breaks,plot=F)$counts 
    entry <- data.frame(tree=rep(t,length(treeCounts)),cds=x,chrom=rep(chr,length(treeCounts)),prop=treeCounts/total)
    cdsFrame <- rbind(cdsFrame,entry)
  }
}

cdsFrame$cds <- round(cdsFrame$cds,2)
cdsFrame$cds<- as.factor(cdsFrame$cds)

#top two
plot <- ggplot(subset(cdsFrame, tree %in% topTwo))+
  geom_boxplot(aes(x=cds,y=prop,col=tree), show.legend=F, outlier.shape=NA) +
  labs(x="Coding Density (ln(Coding Bases per 10 kb Window))",y="Fraction of Trees")+
  scale_x_discrete(labels= c("0-1","1-4","4-5","5-6","6-7","7-8",">8"))+
  #annotate("text",label=nNum,x=seq(1,length(breaks)-1),y=.8)+
  theme(text=element_text(size=20))+
  scale_color_manual(values = colorOrder)
plot

ggsave(plot, filename = paste0(out,"treesByCodingDensity_topTwo.pdf"), width = 10, height=10)

#all trees
plot <- ggplot(subset(cdsFrame, tree %in% topTrees))+
  geom_boxplot(aes(x=cds,y=prop,col=tree), show.legend=F, outlier.shape=NA) +
  labs(x="Coding Density (ln(Coding Bases per 10 kb Window))",y="Fraction of Trees")+
  scale_x_discrete(labels= c("0-1","1-4","4-5","5-6","6-7","7-8",">8"))+
  #annotate("text",label=nNum,x=seq(1,length(breaks)-1),y=.8)+
  theme(text=element_text(size=20))+
  scale_color_manual(values = colorOrder)
plot

ggsave(plot, filename = paste0(out,"treesByCodingDensity_all.pdf"), width = 10, height=10)

for(bin in unique(cdsFrame$cds)){
  print(bin)
  print(t.test(paired=T, x=subset(cdsFrame, cds==bin & tree==topTwo[1])$prop,y=subset(cdsFrame, cds==bin & tree==topTwo[2])$prop))
}





########## control plots Figures S7.8-S7.13 ###########
#recombination rate vs chromosomal position
recombVsPos <- ggplot(data=treeRate) +
  geom_point(aes(x=.5-distFromCenter,y=LD), alpha=.5) +
  geom_smooth(aes(x=.5-distFromCenter,y=LD)) +
  theme(text=element_text(size=20))+
  labs(x="Distance from Edge",y="Recombination Rate (cross-based)")
recombVsPos

#gene density vs chromosomal position
CDSVsPos <- ggplot(data=treeRate) +
  geom_point(aes(x=.5-distFromCenter,y=CDS), alpha=.5) +
  geom_smooth(aes(x=.5-distFromCenter,y=CDS))+
  theme(text=element_text(size=20))+
  labs(x="Distance from Edge",y="Number of Coding Base Pairs per Window")
CDSVsPos

#gene density (per cM) vs chromosomal position
CDSVsPos <- ggplot(data=treeRate) +
  geom_point(aes(x=.5-distFromCenter,y=CDSperCM), alpha=.5) +
  geom_smooth(aes(x=.5-distFromCenter,y=CDSperCM))+
  theme(text=element_text(size=20))+
  labs(x="Distance from Edge",y="Number of Coding Base Pairs per cM")
CDSVsPos

#gene density vs recombination rate
CDSVsRecomb <- ggplot(data=subset(treeRate)) +
  geom_point(aes(y=LD,x=CDS), alpha=.5) +
  geom_smooth(aes(y=LD,x=CDS))+
  theme(text=element_text(size=20))+
  labs(x="Number of Coding Base Pairs per Window",y="Recombination Rate")
CDSVsRecomb

#gene density (CDSbp/CM) vs recombination rate 
CDScMVsRecomb <- ggplot(data=subset(treeRate)) +
  geom_point(aes(y=LD,x=CDSperCM), alpha=.5) +
  geom_smooth(aes(y=LD,x=CDSperCM))+
  theme(text=element_text(size=20))+
  labs(x="Coding Base Pairs per cM",y="Recombination Rate")
CDScMVsRecomb

#log gene density vs recombination rate
logCDScMVsRecomb <- ggplot(data=subset(treeRate)) +
  geom_point(aes(y=LD,x=logCDSdens), alpha=.5) +
  geom_smooth(aes(y=LD,x=logCDSdens))+
  theme(text=element_text(size=20))+
  labs(x="ln(Coding Base Pairs per 10 kb region)",y="Recombination Rate (cM/MB)")
logCDScMVsRecomb

#log gene density vs recombination rate
logCDScMVsRecomb <- ggplot(data=subset(treeRate)) +
  geom_point(aes(y=LD,x=log(CDSperCM+1)), alpha=.5) +
  geom_smooth(aes(y=LD,x=log(CDSperCM+1)))+
  theme(text=element_text(size=20))+
  labs(x="ln(Coding Base Pairs per cM)",y="Recombination Rate (cM/MB")
logCDScMVsRecomb

#gene density vs chromosome size
GVSframe <- data.frame(chrom=HeratoGenome$chromosome, size=HeratoGenome$end)
numCoding <- c()
for(x in GVSframe$chrom){
  numCoding <- c(numCoding,sum(subset(treeRate, chrom==x)$CDS))
}
GVSframe$pctCoding=numCoding/GVSframe$size
CDSvsSize <- ggplot(data=GVSframe) +
  geom_point(aes(x=size/1000000,y=pctCoding)) +
  geom_text_repel(aes(label=chrom,y=0.1,x=size/1000000),direction="x",nudge_y=.02, size=5)+
  theme(text=element_text(size=20))+
  labs(x="Chromosome Size (MB)",y="Percent Coding Bases")
CDSvsSize

treeRate=subset(treeRate, ! chrom %in% c("chr21"))#,"chr15","chr2"))
#topology vs gene density
treeRateSub <- subset(treeRate, tree %in% c( "tree0","tree3","tree2","tree8","tree4","tree6","tree1","tree12"))
treeRateSub$tree <- factor(treeRateSub$tree, levels=c( "tree3","tree0","tree2","tree8","tree1","tree12","tree4","tree6"))
CDSVsTree <- ggplot(data=treeRateSub) +
  geom_boxplot(aes(x=tree, y=CDSperCM,col=tree), notch=T, show.legend=F) +
  #geom_jitter(aes(x=tree, y=CDS,col=tree), alpha=.1)+
  scale_color_manual(values= c(colors[6],colors[2],colors[3],colors[4],colors[8],colors[1],colors[7],colors[5]))+
  scale_x_discrete(labels=c("Tree1", "Tree2", "Tree3", "Tree4", "Tree5", "Tree6", "Tree7", "Tree8"))+
  theme(text=element_text(size=20))+
  labs(y="Number of Coding Base Pairs per Window") 
CDSVsTree

#Zero genes vs topology
ZeroCDSVsTree <- ggplot(data=subset(treeRateSub, CDS==0), aes(tree)) +
  geom_bar(aes(y=..count../sum(..count..),fill=tree)) +
  #geom_jitter(aes(x=tree, y=CDS,col=tree), alpha=.1)+
  scale_fill_manual(values= c(colors[2],colors[6],colors[3],colors[4],colors[7],colors[5],colors[8],colors[1]))+
  labs(y="Percentage of Windows with Zero Coding Bases that Recover Tree",x="Tree")
ZeroCDSVsTree

#nonZero genes vs topology
nonZeroCDSVsTree <- ggplot(data=subset(treeRateSub, CDS!=0), aes(tree)) +
  geom_bar(aes(y=..count../sum(..count..),fill=tree)) +
  #geom_jitter(aes(x=tree, y=CDS,col=tree), alpha=.1)+
  scale_fill_manual(values= c(colors[2],colors[6],colors[3],colors[4],colors[7],colors[5],colors[8],colors[1]))+
  labs(y="Percentage of Windows with Non-Zero Coding Bases that Recover Tree",x="Tree")
nonZeroCDSVsTree



