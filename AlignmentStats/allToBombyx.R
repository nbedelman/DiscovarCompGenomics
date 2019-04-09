library(ggplot2)

exonOverlaps=read.csv("AlignmentStats/alignmentDepth_allGenomes_bombyxRef_exons.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","exonScaf","exonLabel","exonLabel2","exonStart","exonEnd","dot","dash","dot2","ID","overlapLength"))
exonOverlapsTotal=sum(exonOverlaps$overlapLength)

exons=read.csv("AlignmentStats/Bmor_cDNA_exon.gff", sep="\t",header=F, col.names=c("scaf","label","exon","start","end","dot","dash","dot2","ID"))
exons$length=exons$end-exons$start
exonsTotal=sum(exons$length)

####

intronOverlaps=read.csv("AlignmentStats/alignmentDepth_allGenomes_bombyxRef_introns.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intronScaf","intronLabel","intronLabel2","intronStart","intronEnd","dot","dash","dot2","ID","overlapLength"))
intronOverlapsTotal=sum(intronOverlaps$overlapLength)

introns=read.csv("AlignmentStats/Bmor_cDNA_intron.gff", sep="\t",header=F, col.names=c("scaf","label","intron","start","end","dot","dash","dot2","ID"))
introns$length=introns$end-introns$start
intronsTotal=sum(introns$length)

###

intergenicOverlaps=read.csv("AlignmentStats/alignmentDepth_allGenomes_bombyxRef_intergenic.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intergenicScaf","intergenicStart","intergenicEnd","overlapLength"))
intergenicOverlapsTotal=sum(intergenicOverlaps$overlapLength)

intergenic=read.csv("AlignmentStats/Bmor_cDNA_intergenic.bed", sep="\t",header=F, col.names=c("scaf","start","end"))
intergenic$length=intergenic$end-intergenic$start
intergenicTotal=sum(intergenic$length)

###

overallOverlaps=read.csv("AlignmentStats/alignmentDepth_allGenomes_bombyxRef_overall.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intergenicScaf","intergenicStart","intergenicEnd","overlapLength"))
overallOverlapsTotal=sum(overallOverlaps$overlapLength)

overall=read.csv("AlignmentStats/Bmor_fullScaffolds.bed", sep="\t",header=F, col.names=c("scaf","start","end"))
overall$length=overall$end-overall$start
overallTotal=sum(overall$length)

###

pctExonCov=exonOverlapsTotal/exonsTotal
pctIntronsCov=intronOverlapsTotal/intronsTotal
pctintergenicCov=intergenicOverlapsTotal/intergenicTotal
pctOverallCov=overallOverlapsTotal/overallTotal

###
df <- data.frame(labels=c("exons","introns","intergenic", "overall"), 
                 align=c(rep("All Genomes to B. mori",4)), 
                 values=c(pctExonCov,pctIntronsCov,pctintergenicCov,pctOverallCov))

df$labels=factor(df$labels, levels=c("exons","introns","intergenic","overall"))

x <- ggplot(df)+
  geom_bar(aes(x=labels,y=values,fill=align), stat="identity",position = "dodge")+
  geom_text(aes(x=labels,y=values+.01,label=round(values,4))) +
  labs(x="",y="Fraction",title="Fraction Aligned")
x


totalGenomeLength=456867856
totalOverlapLength=exonOverlapsTotal+intronOverlapsTotal+intergenicOverlapsTotal
