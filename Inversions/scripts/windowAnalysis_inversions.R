############## import Libraries ###############

library(devtools)
library("karyoploteR")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plotly)
library(tidyr)
library(knitr)
library(RColorBrewer)

setwd("~/Documents/Mallet_Lab/DiscovarCompGenomics/")

########### Define Data Files #############

ABBABABA="introgression/eratoClade/50KBsliding/ABBABABAsummary.csv"
TREES="introgression/eratoClade/50KBsliding/foundTrees.csv"
DXY="introgression/eratoClade/50KBsliding/distanceMatrix.csv"
out="eratoClade_50KBsliding"

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
Dstat$segment <- as.character(Dstat$segment)
Dstat$chrom <- as.character(unlist(lapply(Dstat$segment,getChrom)))

########  50KB windows ############
Dstat$end <- as.numeric(unlist(lapply(Dstat$segment,getEnd)))*25000
Dstat$start <- Dstat$end-24999
Dstat <- subset(Dstat, !grepl("reduced", segment))
Dstat <- subset(Dstat, !grepl("multi", segment))
Dstat <- Dstat[,c(7,9,8,1,2,3,4,5,6)]

DstatLarge <- subset(Dstat,alnLength>10000 & alnLength<60000 & chrom!="chr0")

Distances <- read.csv(DXY, header=FALSE,col.names = c("segment","EtalHtel","EtalHhsa","EtalHhim","EtalHera","EtalHdem","EtalHsar",
                                                      "HtelHhsa","HtelHhim","HtelHera","HtelHdem","HtelHsar","HhsaHhim","HhsaHera",
                                                      "HhsaHdem","HhsaHsar","HhimHera","HhimHdem","HhimHsar","HeraHdem","HeraHsar",
                                                      "HdemHsar"))
Distances$segment <- as.character(Distances$segment)
Distances$chrom <- as.character(unlist(lapply(Distances$segment,getChrom)))
Distances$end <- as.numeric(unlist(lapply(Distances$segment,getEnd)))*25000
Distances$start <- Distances$end-24999
Distances$overallAverage <- rowMeans(Distances[,2:22])

########### combine data ############

allData <- inner_join(x=DstatLarge, y=trees, by="segment")
allData <- inner_join(x=allData, y=Distances, by="segment")
names <- c("chrom","start","end",colnames(allData)[4:length(colnames(allData))])
colnames(allData) <- names
allData <- subset(allData, alnLength>10000 & alnLength < 60000)

allData$start <- as.character(allData$start)
allData$end <- as.character(allData$end)

ABGranges <- makeGRangesFromDataFrame(allData,keep.extra.columns = TRUE)
ABGranges <- sortSeqlevels(ABGranges)
ABGranges <- sort(ABGranges)


####### Chromosome 15 inversion ##############
tel1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Htel.EI_a_scaffold_13664.chroms.bed")
tel1_15Thick <- toDataframe(tel1_15)
tel1_15Thick <- makeGRangesFromDataFrame(df=tel1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
tel1_15Thick$itemRGB <- tel1_15$itemRgb
tel1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Htel.EI_a_scaffold_13664.chroms.thin.bed")

tel2_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Htel.EI_a_scaffold_33803.chroms.bed")
tel2_15Thick <- toDataframe(tel2_15)
tel2_15Thick <- makeGRangesFromDataFrame(df=tel2_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
tel2_15Thick$itemRGB <- tel2_15$itemRgb
tel2_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Htel.EI_a_scaffold_33803.chroms.thin.bed")

dem1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hdem.EI_a_scaffold_13937.chroms.bed")
dem1_15Thick <- toDataframe(dem1_15)
dem1_15Thick <- makeGRangesFromDataFrame(df=dem1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
dem1_15Thick$itemRGB <- dem1_15$itemRgb
dem1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hdem.EI_a_scaffold_13937.chroms.thin.bed")

dem2_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hdem.EI_a_scaffold_15584.chroms.bed")
dem2_15Thick <- toDataframe(dem2_15)
dem2_15Thick <- makeGRangesFromDataFrame(df=dem2_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
dem2_15Thick$itemRGB <- dem2_15$itemRgb
dem2_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hdem.EI_a_scaffold_15584.chroms.thin.bed")


hsa1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hhsa.EI_a_scaffold_14250.chroms.bed")
hsa1_15Thick <- toDataframe(hsa1_15)
hsa1_15Thick <- makeGRangesFromDataFrame(df=hsa1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
hsa1_15Thick$itemRGB <- hsa1_15$itemRgb
hsa1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hhsa.EI_a_scaffold_14250.chroms.thin.bed")

hsa2_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hhsa.EI_a_scaffold_14429.chroms.bed")
hsa2_15Thick <- toDataframe(hsa2_15)
hsa2_15Thick <- makeGRangesFromDataFrame(df=hsa2_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
hsa2_15Thick$itemRGB <- hsa2_15$itemRgb
hsa2_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hhsa.EI_a_scaffold_14429.chroms.thin.bed")

sar1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hsar.EI_a_scaffold_8872.chroms.bed")
sar1_15Thick <- toDataframe(sar1_15)
sar1_15Thick <- makeGRangesFromDataFrame(df=sar1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
sar1_15Thick$itemRGB <- sar1_15$itemRgb
sar1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hsar.EI_a_scaffold_8872.chroms.thin.bed")

era1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/noInv/HeraRef.chr15_5_noInv.chroms.bed")
era1_15Thick <- toDataframe(era1_15)
era1_15Thick <- makeGRangesFromDataFrame(df=era1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
era1_15Thick$itemRGB <- era1_15$itemRgb
era1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/noInv/HeraRef.chr15_5_noInv.chroms.thin.bed")

him1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/noInv/Hhim.EI_a_scaffold_38868_noInv.chroms.bed")
him1_15Thick <- toDataframe(him1_15)
him1_15Thick <- makeGRangesFromDataFrame(df=him1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
him1_15Thick$itemRGB <- him1_15$itemRgb

Etal1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/noInv/Etal.EI_a_scaffold_91_noInv.chroms.bed")
Etal1_15Thick <- toDataframe(Etal1_15)
Etal1_15Thick <- makeGRangesFromDataFrame(df=Etal1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
Etal1_15Thick$itemRGB <- Etal1_15$itemRgb
Etal1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/noInv/Etal.EI_a_scaffold_91_noInv.chroms.thin.bed")

Etal2_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/noInv/Etal.EI_a_scaffold_41857_noInv.chroms.bed")
Etal2_15Thick <- toDataframe(Etal2_15)
Etal2_15Thick <- makeGRangesFromDataFrame(df=Etal2_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
Etal2_15Thick$itemRGB <- Etal2_15$itemRgb
Etal2_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/noInv/Etal.EI_a_scaffold_41857_noInv.chroms.thin.bed")

num1_15 <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hnum.EI_a_scaffold_16807.chroms_repeatRemoved.bed")
num1_15Thick <- toDataframe(num1_15)
num1_15Thick <- makeGRangesFromDataFrame(df=num1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
num1_15Thick$itemRGB <- num1_15$itemRgb
num1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hnum.EI_a_scaffold_16807.chroms_repeatRemoved.thin.bed")


par1_15<- toGRanges("introgression/chr15_fineScale/relevant_beds/Hpar.EI_a_scaffold_39254.chroms.bed")
par1_15Thick <- toDataframe(par1_15)
par1_15Thick <- makeGRangesFromDataFrame(df=par1_15Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
par1_15Thick$itemRGB <- par1_15$itemRgb
par1_15Thin <- toGRanges("introgression/chr15_fineScale/relevant_beds/Hpar.EI_a_scaffold_39254.chroms.thin.bed")



####### Chromosome 2 inversion ##########
hsa1_2 <- toGRanges("datasets/chr2_fineScale/relevant beds/Hhsa.EI_a_scaffold_11779_singleCopy.goodSegs.bed.condensed.chroms.bed")
hsa1_2Thick <- toDataframe(hsa1_2)
hsa1_2Thick <- makeGRangesFromDataFrame(df=hsa1_2Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
hsa1_2Thick$itemRGB <- hsa1_2$itemRgb
hsa1_2Thin <- toGRanges("datasets/chr2_fineScale/relevant beds/Hhsa.EI_a_scaffold_11779_singleCopy.goodSegs.bed.condensed.chroms.thin.bed")

era5_2 <- toGRanges("datasets/chr2_fineScale/relevant beds/HeraRef.Herato_chr2_5_singleCopy.goodSegs.bed.condensed.chroms.bed")
era5_2Thick <- toDataframe(era5_2)
era5_2Thick <- makeGRangesFromDataFrame(df=era5_2Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
era5_2Thick$itemRGB <- era5_2$itemRgb
era5_2Thin <- toGRanges("datasets/chr2_fineScale/relevant beds/HeraRef.Herato_chr2_5_singleCopy.goodSegs.bed.condensed.chroms.thin.bed")

eraDisco_2 <- toGRanges("datasets/chr2_fineScale/relevant beds/HeraDisco.EI_a_scaffold_22461_singleCopy.bed.condensed.chroms.bed")
eraDisco_2Thick <- toDataframe(eraDisco_2)
eraDisco_2Thick <- makeGRangesFromDataFrame(df=eraDiscoThick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
eraDisco_2Thick$itemRGB <- eraDisco_2$itemRgb
eraDisco_2Thin <- toGRanges("datasets/chr2_fineScale/relevant beds/HeraDisco.EI_a_scaffold_22461_singleCopy.bed.condensed.chroms.thin.bed")


dem1_2 <- toGRanges("datasets/chr2_fineScale/relevant beds/Hdem.EI_a_scaffold_16044_singleCopy.bed.condensed.chroms.bed")
dem1_2Thick <- toDataframe(dem1_2)
dem1_2Thick <- makeGRangesFromDataFrame(df=dem1_2Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
dem1_2Thick$itemRGB <- dem1_2$itemRgb
dem1_2Thin <- toGRanges("datasets/chr2_fineScale/relevant beds/Hdem.EI_a_scaffold_16044_singleCopy.bed.condensed.chroms.thin.bed")

sar1_2 <- toGRanges("datasets/chr2_fineScale/relevant beds/Hsar.EI_a_scaffold_14469_singleCopy.bed.condensed.chroms.bed")
sar1_2Thick <- toDataframe(sar1_2)
sar1_2Thick <- makeGRangesFromDataFrame(df=sar1_2Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
sar1_2Thick$itemRGB <- sar1_2$itemRgb
sar1_2Thin <- toGRanges("datasets/chr2_fineScale/relevant beds/Hsar.EI_a_scaffold_14469_singleCopy.bed.condensed.chroms.thin.bed")

etal1_2 <- toGRanges("datasets/chr2_fineScale/relevant beds/Etal.EI_a_scaffold_46153_singleCopy.bed.condensed.chroms.bed")
etal1_2Thick <- toDataframe(etal1_2)
etal1_2Thick <- makeGRangesFromDataFrame(df=etal1_2Thick,seqnames.field = "chr", start.field = "thick.start",end.field="thick.end")
etal1_2Thick$itemRGB <- etal1_2$itemRgb
etal1_2Thin <- toGRanges("datasets/chr2_fineScale/relevant beds/Etal.EI_a_scaffold_46153_singleCopy.bed.condensed.chroms.thin.bed")



########### Genes ########

genes <- read.delim("data/Hmel2.genes.chroms.bed", header=F, col.names=c("seqnames","start","end","name","score","strand","thickStart","thickEnd","color"))
genes <- makeGRangesFromDataFrame(df=genes, seqnames.field = "seqnames",start.field = "start",end.field = "end", strand.field = "strand", keep.extra.columns = T)
cortex <- genes[which(grepl("HMEL000025",genes$name)),]
cortex$name="Cortex"



####### Plot data to H. melpomene###########
chr15Region=toGRanges(data.frame("chr15", 1.3e6, 2e6))
chr15BP1=toGRanges(data.frame("chr15", 1.3e6, 1.5e6))
chr15BP2=toGRanges(data.frame("chr15", 1.7e6, 1.95e6))
chr15RegionPlus=toGRanges(data.frame("chr15", 1e6, 2.5e6))
chr2Region=toGRanges(data.frame("chr2", 3.6e6, 6.5e6))

pp <- getDefaultPlotParams(plot.type=5)
pp$ideogramheight=0


  #Set up the plot
#to switch between regions, change the chromosome and zoom region here.
  kp <- plotKaryotype(genome=Hmel2pt5GenomeGrange,plot.type=5,cex=2,plot.params = pp,chromosome = c("chr15"), zoom=chr15Region)
  kpAddBaseNumbers(kp,tick.dist = 100000, minor.ticks = T, minor.tick.dist = 10000, cex=2)
  
  ############# Plot Genes ############
  # kpDataBackground(kp, data.panel = 2, r0=0.07, r1=0.1)
  # kpRect(kp,genes,y0=0.3,y1=0.7, r0=0.07, r1=0.1, col="black")
  # kpRect(kp,cortex,y0=0.3,y1=0.7, r0=0.07, r1=0.1, col="gold")


  colors=brewer.pal(8,"Paired")
  ############# Plot trees ############
  #window trees
  kpDataBackground(kp, data.panel = 2,r0=0.01, r1=.05, col="black")
  kpRect(kp, subset(ABGranges,tree=="tree2"),y0=0,y1=1,r0=0.01, r1=.05, col=colors[2], border=NA,data.panel = 2)
  kpRect(kp, subset(ABGranges,tree=="tree1"),y0=0,y1=1,r0=0.01, r1=.05, col=colors[6], border=NA,data.panel = 2)
  kpRect(kp, subset(ABGranges,tree=="tree6"),y0=0,y1=1,r0=0.01, r1=.05, col=colors[3], border=NA,data.panel = 2)
  kpRect(kp, subset(ABGranges,tree=="tree0"),y0=0,y1=1,r0=0.01, r1=.05, col=colors[8], border=NA,data.panel = 2)
  kpRect(kp, subset(ABGranges,tree=="tree3"),y0=0,y1=1,r0=0.01, r1=.05, col=colors[5], border=NA,data.panel = 2)
  kpRect(kp, subset(ABGranges,tree=="tree8"),y0=0,y1=1,r0=0.01, r1=.05, col=colors[7], border=NA,data.panel = 2)
  kpRect(kp, subset(ABGranges,tree=="tree11"),y0=0,y1=1,r0=0.01, r1=.05, col=colors[4], border=NA,data.panel = 2)
  kpRect(kp, subset(ABGranges,! tree %in% c("tree2","tree1","tree6","tree0","tree3","tree8","tree11")),
         y0=0,y1=1,r0=0.01, r1=.05, col="grey", border=NA,data.panel = 2)
  
  # # ################## Chr15_inversion #############
  ########## controls ############# 
  #erato
  kpDataBackground(kp, r0=0.12, r1=0.18)
  kpArrows(kp,era1_15Thin,r0=0.12, r1=0.18, y0=.5,y1=.5,length=0,lty=2, lwd=2)
  kpRect(kp, data=era1_15Thick, y0=0.25, y1=0.74,r0=0.12, r1=0.18,col="black", clipping = T)
  
  #eueides
  kpDataBackground(kp, r0=0.20, r1=0.32)
  kpArrows(kp,Etal2_15Thin,r0=0.20, r1=0.32, y0=.25,y1=.25,length=0,lty=2, lwd=2)
  kpRect(kp, data=Etal2_15Thick, y0=0.1, y1=.4,r0=0.20, r1=0.32,col="black")
  
  kpArrows(kp,Etal1_15Thin,r0=0.20, r1=0.32, y0=.75,y1=.75,length=0,lty=2, lwd=2)
  kpRect(kp, data=Etal1_15Thick, y0=0.6, y1=.9,r0=0.20, r1=0.32,col="black")
  
  ########### inversions ###########
  
  #telesiphe
  kpDataBackground(kp, r0=0.34, r1=0.46)
  kpArrows(kp,tel1_15Thin,r0=0.34, r1=0.46, y0=.25,y1=.25,length=0,lty=2, lwd=2)
  kpRect(kp, data=tel1_15Thick, y0=0.1, y1=.4,r0=0.34, r1=0.46,col="black")

  kpArrows(kp,tel2_15Thin,r0=0.34, r1=0.46, y0=.75,y1=.75,length=0,lty=2, lwd=2)
  kpRect(kp, data=tel2_15Thick,  y0=0.6, y1=.9,r0=0.34, r1=0.46, col="black")

  #hecalesia
  kpDataBackground(kp, r0=0.48, r1=0.60)
  kpArrows(kp,hsa1_15Thin,r0=0.48, r1=0.60, y0=.25,y1=.25,length=0,lty=2, lwd=2)
  kpRect(kp, data=hsa1_15Thick, y0=0.1, y1=.4,r0=0.48, r1=0.60,col="black")
  
  kpArrows(kp,hsa2_15Thin,r0=0.48, r1=0.60, y0=.75,y1=.75,length=0,lty=2, lwd=2)
  kpRect(kp, data=hsa2_15Thick,  y0=0.6, y1=.9,r0=0.48, r1=0.60, col="black")
  
  #demeter
  kpDataBackground(kp, r0=0.62, r1=.74)
  kpArrows(kp,dem1_15Thin,r0=0.62, r1=.74, y0=.25,y1=.25,length=0,lty=2, lwd=2)
  kpRect(kp, data=dem1_15Thick, y0=0.1, y1=.4,r0=0.62, r1=.74,col="black")
  
  kpArrows(kp,dem2_15Thin,r0=0.62, r1=.74, y0=.75,y1=.75,length=0,lty=2, lwd=2)
  kpRect(kp, data=dem2_15Thick,  y0=0.6, y1=.9,r0=0.62, r1=.74, col="black")

  #sara
  kpDataBackground(kp, r0=0.76, r1=.82)
  kpArrows(kp,sar1_15Thin,r0=0.76, r1=.82, y0=.5,y1=.5,length=0,lty=2, lwd=2)
  kpRect(kp, data=sar1_15Thick, y0=0.25, y1=.75,r0=0.76, r1=.82,col="black")

  
  
  #numata
  kpDataBackground(kp, r0=0.84, r1=.90)
  kpArrows(kp,num1_15Thin,r0=0.84, r1=.90, y0=.5,y1=.5,length=0,lty=2, lwd=2)
  kpRect(kp, data=num1_15Thick, y0=0.25, y1=.75,r0=0.84, r1=.90,col="black")

  
  # # ################## Chr2_inversion #############
  ### controls
  
  kpDataBackground(kp, r0=0.12, r1=0.18)
  kpArrows(kp,etal1Thin,r0=0.12, r1=0.18, y0=.5,y1=.5,length=0,lty=2, lwd=2)
  kpRect(kp, data=etal1Thick, y0=0.25, y1=0.74,r0=0.12, r1=0.18,col="black", clipping = T)
  
  kpDataBackground(kp, r0=0.20, r1=0.26)
  kpArrows(kp,sar1Thin,r0=0.20, r1=0.26, y0=.5,y1=.5,length=0,lty=2, lwd=2)
  kpRect(kp, data=sar1Thick, y0=0.25, y1=0.74,r0=0.20, r1=0.26,col="black", clipping = T)
  
  ### inversions
  
  kpDataBackground(kp, r0=0.28, r1=0.34)
  kpArrows(kp,era5Thin,r0=0.28, r1=0.34,  y0=.5,y1=.5,length=0,lty=2, lwd=2)
  kpRect(kp, data=era5Thick, y0=0.25, y1=0.74,r0=0.28, r1=0.34,col="black", clipping = T)

  kpDataBackground(kp, r0=0.36, r1=0.42)
  kpArrows(kp,hsa1Thin,r0=0.36, r1=0.42, y0=.5,y1=.5,length=0,lty=2, lwd=2)
  kpRect(kp, data=hsa1Thick, y0=0.25, y1=0.74,r0=0.36, r1=0.42,col="black", clipping = T)

