############## import Libraries ###############

#source("https://bioconductor.org/biocLite.R")
#install.packages("knitr")
library(devtools)
#install_github("bernatgel/karyoploteR")
# biocLite("karyoploteR")
library("karyoploteR")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
# library(plotly)
# library(tidyr)
# library(plot3Drgl)
# library(plot3D)
# library(rgl)
# library(knitr)
library(RColorBrewer)

########### Define Data Files #############

ABBABABA="introgression/melpomeneClade/ABBABABAsummary.csv"
TREES="introgression/melpomeneClade/foundTrees.csv"
DXY="introgression/eratoClade/50KBsliding/distanceMatrix.csv"
out="melpomeneClade_50KBSliding"

# for (i in 1:length(args)) {
#   eval (parse (text = args[[i]] ))
# }

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


################Melpomene genome setup################
Hmel2pt5 <- read.delim("introgression/orderChroms/Hmel2.transitions.tsv") 
Hmel2pt5$Chromosome <- paste0("chr",Hmel2pt5$Chromosome)
Hmel2pt5Genome <- read.delim("introgression/orderChroms/Hmel2.transitions_fullChroms.bed", col.names=c("seqnames","start","end","chromosome"))
Hmel2pt5Genome$chr <-  unlist(strsplit(as.character(Hmel2pt5Genome$chromosome),"chr"))[c(FALSE,TRUE)]
Hmel2pt5Cyto <- data.frame(seqnames=Hmel2pt5$Chromosome, Hmel2pt5$ChromStart,Hmel2pt5$ChromEnd, strand=Hmel2pt5$Orientation, label=Hmel2pt5$Hmel2Scaffold)

Hmel2pt5Genome <- Hmel2pt5Genome[order(as.numeric(Hmel2pt5Genome[,5])),]

Hmel2pt5GenomeGrange <- makeGRangesFromDataFrame(Hmel2pt5Genome, seqnames.field = "seqnames",keep.extra.columns = TRUE)
Hmel2pt5CytoGrange <- makeGRangesFromDataFrame(Hmel2pt5Cyto)
trees <- read.csv(TREES, header=FALSE, col.names = c("segment","tree"))
trees$segment <- as.character(trees$segment)

Dstat <- read.csv(ABBABABA, header=FALSE,col.names = c("segment","alnLength","snps","biallelic","ABBA","BABA"))
#Dstat <- read.csv(ABBABABA, header=FALSE,col.names = c("segment","alnLength","snps","biallelic"))
Dstat$segment <- as.character(Dstat$segment)
Dstat$chrom <- as.character(unlist(lapply(Dstat$segment,getChrom)))

########  50KB windows ############
# 
Dstat$end <- as.numeric(unlist(lapply(Dstat$segment,getEnd)))*25000
Dstat$start <- Dstat$end-24999
Dstat <- subset(Dstat, !grepl("reduced", segment))
Dstat <- subset(Dstat, !grepl("multi", segment))
Dstat <- Dstat[,c(7,9,8,1,2,3,4,5,6)]
# DstatGranges <- makeGRangesFromDataFrame(df=Dstat,keep.extra.columns = TRUE)
# DstatGranges <- sortSeqlevels(DstatGranges)
# DstatGranges <- sort(DstatGranges)
DstatLarge <- subset(Dstat,alnLength>5000 & alnLength<60000 & chrom!="chr0")

Distances <- read.csv(DXY, header=FALSE,col.names = c("segment","EtalHtel","EtalHhsa","EtalHhim","EtalHera","EtalHdem","EtalHsar",
                                                      "HtelHhsa","HtelHhim","HtelHera","HtelHdem","HtelHsar","HhsaHhim","HhsaHera",
                                                      "HhsaHdem","HhsaHsar","HhimHera","HhimHdem","HhimHsar","HeraHdem","HeraHsar",
                                                      "HdemHsar"))
Distances$segment <- as.character(Distances$segment)
Distances$chrom <- as.character(unlist(lapply(Distances$segment,getChrom)))
Distances$end <- as.numeric(unlist(lapply(Distances$segment,getEnd)))*25000
Distances$start <- Distances$end-24999
Distances$overallAverage <- rowMeans(Distances[,2:22])
# Distances$alnLength <- as.numeric(unlist(lapply(Distances$segment,getAlnLength)))
# Distances=subset(Distances, alnLength>10000& alnLength<60000)
########### combine data ############
# allData <- trees
# allData <- Dstat
allData <- inner_join(x=Dstat, y=trees, by="segment")
# allData <- inner_join(x=allData, y=Distances, by="segment")
# names <- c("chrom","start","end",colnames(allData)[4:length(colnames(allData))])
# colnames(allData) <- names
# allData <- subset(allData, alnLength>5000 & alnLength < 60000)

##10KB every 50 trial
# trees$chrom <- as.character(unlist(lapply(trees$segment,getChrom)))
# trees$end <- as.numeric(unlist(lapply(trees$segment,getEnd)))*50000
# trees$start <- trees$end-49999



###############
allData$start <- as.character(allData$start)
allData$end <- as.character(allData$end)

bestTrees <- names(summary(allData$tree))[1:50]
allData <- subset(allData, tree %in% bestTrees)
allData$tree <- as.character(allData$tree)
allData$tree <- factor(allData$tree, levels=bestTrees)

ABGranges <- makeGRangesFromDataFrame(allData,keep.extra.columns = TRUE)
ABGranges <- sortSeqlevels(ABGranges)
ABGranges <- sort(ABGranges)

####### Plot data to H. melpomene###########
chr15Region=toGRanges(data.frame("chr15", 1.3e6, 2e6))
chr15BP1=toGRanges(data.frame("chr15", 1.3e6, 1.5e6))
chr15BP2=toGRanges(data.frame("chr15", 1.7e6, 1.95e6))
chr15RegionPlus=toGRanges(data.frame("chr15", 1e6, 2.5e6))
chr2Region=toGRanges(data.frame("chr2", 4e6, 6.5e6))

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight=0


for(chr in c(as.character(Hmel2pt5Genome$seqnames),"")){
  
tiff(filename = paste0(out,"eratoClade_trees_toHmel.tiff"), width=2000,height=1000)
#Set up the plot
kp <- plotKaryotype(genome=Hmel2pt5GenomeGrange,plot.type=1,cex=2,plot.params = pp)#,chromosome = c("chr9"))#, zoom=chr15Region)#cytobands = Hmel2pt5CytoGrange, cex=2)#, chromosome = c("chr15"))#, zoom = chr15Region)
#kp <- plotKaryotype(genome=HeratoGenomeGrange,plot.type=1,cex=2,plot.params = pp)#,chromosome = c("chr2"))#, zoom=chr15Region)#cytobands = Hmel2pt5CytoGrange, cex=2)#, chromosome = c("chr15"))#, zoom = chr15Region)
  
#kpAddBaseNumbers(kp,tick.dist = 1000000, minor.ticks = T, minor.tick.dist = 100000, cex=2)

colors=brewer.pal(8,"Paired")
getPalette <-  colorRampPalette(brewer.pal(9, "Set1"))
thisPal <- getPalette(75)
palette(thisPal)

############# Plot trees ############
#window trees
kpDataBackground(kp, data.panel = 1,r0=0, r1=.75, col="black")
kpRect(kp, data=ABGranges,y0=0,y1=1,r0=0, r1=.75, col=ABGranges$tree, border=NA)
legend("topright", legend=bestTrees,
       fill=thisPal, pch=15, cex=.80)
kpRect(kp, subset(ABGranges,tree=="tree44"),y0=0,y1=1,r0=0, r1=.75, col=thisPal[1], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree11"),y0=0,y1=1,r0=0,r1=.75, col=colors[6], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree31"),y0=0,y1=1,r0=0, r1=.75, col=colors[3], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree35"),y0=0,y1=1,r0=0, r1=.75, col=colors[8], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree26"),y0=0,y1=1,r0=0, r1=.75, col=colors[5], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree47"),y0=0,y1=1,r0=0, r1=.75, col=colors[7], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree91"),y0=0,y1=1,r0=0, r1=.75, col=colors[4], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree92"),y0=0,y1=1,r0=0.03, r1=0.75, col=colors[1], border=NA)
#kpRect(kp, subset(ABGranges,! tree %in% c("tree2","tree1","tree6","tree0","tree3","tree8","tree11")),
#       y0=0,y1=1,r0=0, r1=.25, col="grey", border=NA)

}

dev.off()

data <- iris
plot(data$Sepal.Length, data$Sepal.Width, col=data$Species)

kpDataBackground(kp, data.panel = 1,r0=0, r1=.75, col="black")
kpRect(kp, subset(ABGranges,tree=="tree0"),y0=0,y1=1,r0=0, r1=.75, col=colors[2], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree1"),y0=0,y1=1,r0=0,r1=.75, col=colors[6], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree2"),y0=0,y1=1,r0=0, r1=.75, col=colors[3], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree8"),y0=0,y1=1,r0=0, r1=.75, col=colors[8], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree9"),y0=0,y1=1,r0=0, r1=.75, col=colors[5], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree6"),y0=0,y1=1,r0=0, r1=.75, col=colors[7], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree7"),y0=0,y1=1,r0=0, r1=.75, col=colors[4], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree12"),y0=0,y1=1,r0=0, r1=0.75, col=colors[1], border=NA)
kpRect(kp, subset(ABGranges,! tree %in% c("tree4","tree5","tree1","tree8","tree9","tree6","tree7","tree12")),
       y0=0,y1=1,r0=0, r1=.75, col="grey", border=NA)



