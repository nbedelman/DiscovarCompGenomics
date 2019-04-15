#source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")
library(devtools)
install_github("bernatgel/karyoploteR")
#biocLite("karyoploteR")
library("karyoploteR")
library(ggplot2)

args=commandArgs(TRUE)

#for (i in 1:length(args)) {
#  eval (parse (text = args[[i]] ))
#}

ABBABABA=args[[1]]
TREES=args[[2]]
out=args[[3]]
mapFile=args[[4]]
fullChroms=args[[5]]

getChrom <- function (x) {
  chrom <- strsplit(x,"_")[[1]][1]
  return(chrom)
}

getEnd <- function (x) {
  end <- strsplit(x,"_")[[1]][2]
  return(end)
}

Hmel2pt5 <- read.delim(mapFile)
Hmel2pt5$Chromosome <- paste0("chr",Hmel2pt5$Chromosome)
Hmel2pt5Genome <- read.delim(fullChroms, col.names=c("seqnames","start","end","chromosome"))
Hmel2pt5Genome$chr <-  unlist(strsplit(as.character(Hmel2pt5Genome$chromosome),"chr"))[c(FALSE,TRUE)]
Hmel2pt5Cyto <- data.frame(seqnames=Hmel2pt5$Chromosome, Hmel2pt5$ChromStart,Hmel2pt5$ChromEnd, strand=Hmel2pt5$Orientation, label=Hmel2pt5$Hmel2Scaffold)

Hmel2pt5Genome <- Hmel2pt5Genome[order(as.numeric(Hmel2pt5Genome[,5])),]

Hmel2pt5GenomeGrange <- makeGRangesFromDataFrame(Hmel2pt5Genome, seqnames.field = "seqnames",keep.extra.columns = TRUE)
Hmel2pt5CytoGrange <- makeGRangesFromDataFrame(Hmel2pt5Cyto)

trees <- read.csv(TREES, header=FALSE, col.names = c("region","tree","overallTree"))
trees$region <- as.character(trees$region)

Dstat <- read.csv(ABBABABA, header=FALSE,col.names = c("segment","alnLength","snps","biallelic","ABBA","BABA"))
Dstat$segment <- as.character(Dstat$segment)
Dstat$chrom <- as.character(unlist(lapply(Dstat$segment,getChrom)))
Dstat$end <- as.numeric(unlist(lapply(Dstat$segment,getEnd)))*25000+25000
Dstat$start <- Dstat$end-49999
Dstat <- subset(Dstat, !grepl("reduced", segment))
Dstat <- subset(Dstat, !grepl("multi", segment))
Dstat$D <- (Dstat$ABBA-Dstat$BABA)/(Dstat$ABBA+Dstat$BABA)

allData <- cbind(Dstat,trees)

DstatLarge <- subset(allData,alnLength>10000 & alnLength<60000 & chrom!="chr0" & ABBA+BABA>20)
#DstatLarge <- subset(Dstat,alnLength>10000 & alnLength<60000 & chrom!="chr0" & ABBA+BABA>20)


grangeInfo <- data.frame(seqNames=DstatLarge$chrom, start=DstatLarge$start, end=DstatLarge$end, value=DstatLarge$D, tree=DstatLarge$tree, ABBA=DstatLarge$ABBA, BABA=DstatLarge$BABA)
#grangeInfo <- data.frame(seqNames=DstatLarge$chrom, start=DstatLarge$start, end=DstatLarge$end, value=DstatLarge$D, ABBA=DstatLarge$ABBA, BABA=DstatLarge$BABA)
ABGranges <- makeGRangesFromDataFrame(df=grangeInfo,keep.extra.columns = TRUE)


jpeg(filename = paste0(out,"_allAlns.jpeg"), width=2000,height=1000)
alignPlot <- ggplot(data=Dstat)+
  geom_histogram(aes(x=alnLength),bins=100)
alignPlot
dev.off()

jpeg(filename = paste0(out,"_largAlns.jpeg"), width=2000,height=1000)
alignPlot <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=alnLength),bins=100)
alignPlot
dev.off()

jpeg(filename = paste0(out,"_SNPs.jpeg"), width=2000,height=1000)
snpPlot <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=snps),bins=100)
snpPlot
dev.off()

jpeg(filename = paste0(out,"_biallelic.jpeg"), width=2000,height=1000)
biallelicSnpPlot <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=biallelic),bins=100)
biallelicSnpPlot
dev.off()

jpeg(filename = paste0(out,"_ABBAplusBABA.jpeg"), width=2000,height=1000)
ABBAplusBABA <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=ABBA+BABA),binwidth=2)
ABBAplusBABA
dev.off()

jpeg(filename = paste0(out,"_ABBAandBABA.jpeg"), width=2000,height=1000)
ABBAandBABA <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=ABBA), binwidth=2, fill="orange", alpha=0.65)+
  geom_histogram(aes(x=BABA), binwidth=2, fill="purple", alpha=0.65)+
  labs(x="ABBA vs BABA")
ABBAandBABA
dev.off()

#ABBAandBABAauto<- ggplot(data=subset(DstatLarge, !grepl("chr21", segment)))+
#  geom_histogram(aes(x=ABBA), binwidth=2, fill="orange", alpha=0.65)+
#  geom_histogram(aes(x=BABA), binwidth=2, fill="purple", alpha=0.65)+
#  labs(x="ABBA vs BABA")
#ABBAandBABAauto

getMeanSE <- function(seqName){
  vec <- subset(grangeInfo, seqNames==seqName)$value
  return(mean_se(vec))
}

x <- as.character(unique(grangeInfo$seqNames))
summaryStats <- data.frame(seqNames=character(),start=numeric(),end=numeric(),Dmin=numeric(), Dmax=numeric(), ABBAmin=numeric(),
ABBAmax=numeric(), BABAmin=numeric(), BABAmax=numeric())
for (chr in x){
  length=Hmel2pt5Genome[which(Hmel2pt5Genome$chromosome==chr),3]
  #DmeanSE=mean_se(subset(grangeInfo, seqNames==chr)$value)
  #AmeanSE=mean_se(subset(grangeInfo, seqNames==chr)$ABBA)
  #BmeanSE=mean_se(subset(grangeInfo, seqNames==chr)$BABA)
  #newVals <- data.frame(seqNames=chr, start=(length/2)-1000000, end=(length/2)+1000000, Dmin=DmeanSE$ymin, Dmax=DmeanSE$ymax,
  #                      ABBAmin=AmeanSE$ymin, ABBAmax=AmeanSE$ymax,BABAmin=BmeanSE$ymin, BABAmax=BmeanSE$ymax)
  Dmean=mean(subset(grangeInfo, seqNames==chr)$value)
  DSD=sd(subset(grangeInfo, seqNames==chr)$value)
  Amean=mean(subset(grangeInfo, seqNames==chr)$ABBA)
  ASD=sd(subset(grangeInfo, seqNames==chr)$ABBA)
  Bmean=mean(subset(grangeInfo, seqNames==chr)$BABA)
  BSD=sd(subset(grangeInfo, seqNames==chr)$BABA)
  newVals <- data.frame(seqNames=chr, start=(length/2)-1000000, end=(length/2)+1000000, Dmin=Dmean-DSD, Dmax=Dmean+DSD,
                        ABBAmin=Amean-ASD, ABBAmax=Amean+ASD,BABAmin=Bmean-BSD, BABAmax=Bmean+BSD)
  summaryStats <- rbind(summaryStats,newVals)
}
summaryGranges <-  makeGRangesFromDataFrame(summaryStats,keep.extra.columns = TRUE)


totalMeanSE <- mean_se(DstatLarge$D)

####### Plot data ###########

for(chr in c(as.character(unique(grangeInfo$seqNames)),"")){

  jpeg(filename = paste0(out,chr,".jpeg"), width=2000,height=1000)
  #Set up the plot
  kp <- plotKaryotype(genome=Hmel2pt5GenomeGrange,plot.type=4, chromosome = c(chr))

  #Plot all data as scatterplot
  dataMin=round(min(ABGranges$value)-0.1,1)
  dataMax=round(max(ABGranges$value)+0.1,1)
  kpDataBackground(kp, data.panel = 1, r0=0, r1=0.35)
  kpAxis(kp, ymin=dataMin, ymax=dataMax, r0=0, r1=0.35, col="gray50", cex=1, numticks = length(seq(dataMin,dataMax,0.1)))
  kpAbline(kp, h=seq(dataMin,dataMax,0.1), col="grey50",ymin=dataMin, ymax=dataMax, r0=0, r1=0.35)
  kpAbline(kp, h=0, col="red",ymin=dataMin, ymax=dataMax, r0=0, r1=0.35, lwd=2)
  kpPoints(kp, data=ABGranges,cex=.5,r0=0, r1=0.35,ymin=dataMin, ymax=dataMax)

  #Plot mean+/- SE for each chromosome
  sumMin=round(min(summaryGranges$Dmin)-0.1,1)
  sumMax=round(max(summaryGranges$Dmax)+0.1,1)
  kpDataBackground(kp, data.panel = 1, r0=0.45, r1=0.6)
  kpAxis(kp, ymin=sumMin, ymax=sumMax, r0=0.45, r1=0.6, col="gray50", cex=1, numticks = length(seq(sumMin,sumMax,0.1)))
  kpAbline(kp, h=seq(sumMin,sumMax,0.1), col="gray50",ymin=sumMin, ymax=sumMax, r0=0.45, r1=0.6)
  kpAbline(kp, h=0, col="red",ymin=sumMin, ymax=sumMax, r0=0.45, r1=0.6,lwd=2)
  kpRect(kp, data=summaryGranges,y0=(summaryGranges$Dmin+summaryGranges$Dmax)/2,y1=summaryGranges$Dmax,r0=0.45, r1=0.6,ymin=sumMin, ymax=sumMax)
  kpRect(kp, data=summaryGranges,y1=(summaryGranges$Dmin+summaryGranges$Dmax)/2,y0=summaryGranges$Dmin,r0=0.45, r1=0.6,ymin=sumMin, ymax=sumMax)
  kpAbline(kp, h=totalMeanSE$y,r0=0.45, r1=0.6,ymin=sumMin, ymax=sumMax, col="blue", lty="dashed")


  #Plot mean +/- SE for numABBA and numBABA for each chromosome
  ABmin=round(min(summaryGranges$ABBAmin, summaryGranges$BABAmin)-10,-1)
  ABmax=round(max(summaryGranges$ABBAmax, summaryGranges$BABAax)+10,-1)
  kpDataBackground(kp, data.panel = 1, r0=0.65, r1=0.95)
  kpAxis(kp, ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col="gray50",cex=1, numticks = length(seq(ABmin,ABmax,10)))
  kpAbline(kp, h=seq(ABmin,ABmax,10), col="gray50",ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95)
  kpRect(kp, data=summaryGranges,y0=(summaryGranges$ABBAmin+summaryGranges$ABBAmax)/2,y1=summaryGranges$ABBAmax,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col=rgb(1,140/255,0, alpha=0.5))
  kpRect(kp, data=summaryGranges,y1=(summaryGranges$ABBAmin+summaryGranges$ABBAmax)/2,y0=summaryGranges$ABBAmin,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col=rgb(1,140/255,0, alpha=0.5))
  kpRect(kp, data=summaryGranges,y0=(summaryGranges$BABAmin+summaryGranges$BABAmax)/2,y1=summaryGranges$BABAmax,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col=rgb(138/255,43/255,226/255, alpha=0.5))
  kpRect(kp, data=summaryGranges,y1=(summaryGranges$BABAmin+summaryGranges$BABAmax)/2,y0=summaryGranges$BABAmin,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col=rgb(138/255,43/255,226/255, alpha=0.5))


  #Plot trees
  kpDataBackground(kp, data.panel = 1, r0=.36, r1=.40, col="black")
  kpRect(kp, subset(ABGranges,tree=="Expected"),y0=0,y1=1,r0=.36, r1=.40, col="skyblue", border=NA)
  kpRect(kp, subset(ABGranges,tree=="ABBA"),y0=0,y1=1,r0=.36, r1=.40, col="orange", border=NA)
  kpRect(kp, subset(ABGranges,tree=="BABA"),y0=0,y1=1,r0=.36, r1=.40, col="purple", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Eueides+Erato"),y0=0,y1=1,r0=.36, r1=.40, col="orange", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Eueides+Melpomene"),y0=0,y1=1,r0=.36, r1=.40, col="purple", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Silvaniform+Cyd/Tim"),y0=0,y1=1,r0=.36, r1=.40, col="orange", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Silvaniform+Melpomene"),y0=0,y1=1,r0=.36, r1=.40, col="purple", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Erato+Melpomene"),y0=0,y1=1,r0=.36, r1=.40, col="orange", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Sara+Melpomene"),y0=0,y1=1,r0=.36, r1=.40, col="purple", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Melpomene+Erato"),y0=0,y1=1,r0=.36, r1=.40, col="orange", border=NA)
  # kpRect(kp, subset(ABGranges,tree=="Cydno+Erato"),y0=0,y1=1,r0=.36, r1=.40, col="purple", border=NA)
  #
  dev.off()

}
