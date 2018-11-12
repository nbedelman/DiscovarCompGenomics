#Script to generate Figure 2A and Supplementary Figure S8.2

############## Import Libraries ###############
library(devtools)
library("karyoploteR")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plotly)
library(tidyr)
library(knitr)

setwd("HeliconiusComparativeGenomics/")

########### Define Data Files #############
#50KB sliding windows - Figure 2A
ABBABABA="datasets/50KBsliding_toErato/ABBABABAsummary.csv"
TREES="datasets/50KBsliding_toErato/allTrees.found.tsv"
DXY="datasets/50KBsliding_toErato/distanceMatrix.csv"

#10KB abutting windows - Figure S8.2
ABBABABA="datasets/10KBAbutting_toErato/ABBABABAsummary.csv"
TREES="datasets/10KBAbutting_toErato/eraClade.allTrees.found.tsv"
DXY="datasets/10KBAbutting_toErato/eraClade.distances.csv"


HeratoGenome <- read.delim("datasets/Herato.transitions_fullChroms.bed", col.names=c("seqnames","start","end","chromosome"), header=F)
out=""

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


################ Erato genome setup ###########
HeratoGenome$chr <-  unlist(strsplit(as.character(HeratoGenome$chromosome),"chr"))[c(FALSE,TRUE)]
HeratoGenome <- HeratoGenome[order(as.numeric(HeratoGenome[,5])),]
HeratoGenomeGrange <- makeGRangesFromDataFrame(HeratoGenome, seqnames.field = "seqnames",keep.extra.columns = TRUE)

############## Read in data #########
trees <- read.csv(TREES, header=FALSE, col.names = c("segment","tree"))
trees$segment <- as.character(trees$segment)


trees$chrom <- as.character(unlist(lapply(trees$segment,getChrom)))

#10KB windows
trees$end <- as.numeric(unlist(lapply(trees$segment,getEnd)))*10000
trees$start <- trees$end-9999

#50KB windos
trees$end <- as.numeric(unlist(lapply(trees$segment,getEnd)))*25000
trees$start <- trees$end-24999

allData <- trees

###############
allData$start <- as.character(allData$start)
allData$end <- as.character(allData$end)

ABGranges <- makeGRangesFromDataFrame(allData,keep.extra.columns = TRUE)
ABGranges <- sortSeqlevels(ABGranges)
ABGranges <- sort(ABGranges)

####### Plot data to H. erato ###########
pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight=0

kp <- plotKaryotype(genome=HeratoGenomeGrange,plot.type=1,cex=2,plot.params = pp)

colors=brewer.pal(8,"Paired")
#Trees are assigned arbitrary numbers in the filtering process, so need to assign the proper colors separately for each dataset
############# Plot trees for 10KB blocks ############
#window trees
kpDataBackground(kp, data.panel = 1,r0=0, r1=.75, col="black")
kpRect(kp, subset(ABGranges,tree=="tree0"),y0=0,y1=1,r0=0, r1=.75, col=colors[2], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree3"),y0=0,y1=1,r0=0,r1=.75, col=colors[6], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree2"),y0=0,y1=1,r0=0, r1=.75, col=colors[3], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree8"),y0=0,y1=1,r0=0, r1=.75, col=colors[4], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree1"),y0=0,y1=1,r0=0, r1=.75, col=colors[7], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree4"),y0=0,y1=1,r0=0, r1=.75, col=colors[8], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree6"),y0=0,y1=1,r0=0, r1=.75, col=colors[1], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree12"),y0=0,y1=1,r0=0, r1=.75, col=colors[5], border=NA)
kpRect(kp, subset(ABGranges,! tree %in% c("tree0","tree3","tree2","tree8","tree4","tree6","tree1","tree12")),
      y0=0,y1=1,r0=0, r1=0.75, col="grey", border=NA)

############# Plot trees for 50KB blocks ############
#window trees
kpDataBackground(kp, data.panel = 1,r0=0, r1=.75, col="black")
kpRect(kp, subset(ABGranges,tree=="tree4"),y0=0,y1=1,r0=0,r1=.75, col=colors[6], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree2"),y0=0,y1=1,r0=0, r1=.75, col=colors[2], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree0"),y0=0,y1=1,r0=0, r1=.75, col=colors[3], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree5"),y0=0,y1=1,r0=0, r1=.75, col=colors[4], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree6"),y0=0,y1=1,r0=0, r1=.75, col=colors[8], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree9"),y0=0,y1=1,r0=0, r1=.75, col=colors[1], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree3"),y0=0,y1=1,r0=0, r1=.75, col=colors[7], border=NA)
kpRect(kp, subset(ABGranges,tree=="tree7"),y0=0,y1=1,r0=0, r1=.75, col=colors[5], border=NA)
kpRect(kp, subset(ABGranges,! tree %in% c("tree4","tree2","tree0","tree5","tree6","tree9","tree3","tree7")),
       y0=0,y1=1,r0=0, r1=0.75, col="grey", border=NA)


