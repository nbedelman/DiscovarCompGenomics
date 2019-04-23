library(ggplot2)

exonOverlaps=read.csv("referenceAlignmentQuality/coverageDepth_melpomene_discoToRef_exonOverlap.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","exonScaf","exonLabel","exonLabel2","exonStart","exonEnd","dot","dash","dot2","ID","overlapLength"))
exonOverlapsTotal=sum(exonOverlaps$overlapLength)

exons=read.csv("referenceAlignmentQuality/Hmel2_exon.gff", sep="\t",header=F, col.names=c("scaf","label","exon","start","end","dot","dash","dot2","ID"))
exons$length=exons$end-exons$start
exonsTotal=sum(exons$length)

####

intronOverlaps=read.csv("referenceAlignmentQuality/coverageDepth_melpomene_discoToRef_intronOverlap.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intronScaf","intronLabel","intronLabel2","intronStart","intronEnd","dot","dash","dot2","ID","overlapLength"))
intronOverlapsTotal=sum(intronOverlaps$overlapLength)

introns=read.csv("referenceAlignmentQuality/Hmel2_introns.gff", sep="\t",header=F, col.names=c("scaf","label","intron","start","end","dot","dash","dot2","ID"))
introns$length=introns$end-introns$start
intronsTotal=sum(introns$length)

###

intergenicOverlaps=read.csv("referenceAlignmentQuality/coverageDepth_melpomene_discoToRef_intergenicOverlap.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intergenicScaf","intergenicStart","intergenicEnd","overlapLength"))
intergenicOverlapsTotal=sum(intergenicOverlaps$overlapLength)

intergenic=read.csv("referenceAlignmentQuality/Hmel2_intergenic.gff", sep="\t",header=F, col.names=c("scaf","start","end"))
intergenic$length=intergenic$end-intergenic$start
intergenicTotal=sum(intergenic$length)

###

pctExonCov=exonOverlapsTotal/exonsTotal
pctIntronsCov=intronOverlapsTotal/intronsTotal
pctintergenicCov=intergenicOverlapsTotal/intergenicTotal


########################## EratoToMelpomene ########################
EMexonOverlaps=read.csv("AlignmentStats/alignmentDepth_eratoRef_to_HmelRef_exonOverlap.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","exonScaf","exonLabel","exonLabel2","exonStart","exonEnd","dot","dash","dot2","ID","overlapLength"))
EMexonOverlapsTotal=sum(EMexonOverlaps$overlapLength)

####

EMintronOverlaps=read.csv("AlignmentStats/alignmentDepth_eratoRef_to_HmelRef_intronOverlap.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intronScaf","intronLabel","intronLabel2","intronStart","intronEnd","dot","dash","dot2","ID","overlapLength"))
EMintronOverlapsTotal=sum(EMintronOverlaps$overlapLength)
###

EMintergenicOverlaps=read.csv("AlignmentStats/alignmentDepth_eratoRef_to_HmelRef_intergenicOverlap.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intergenicScaf","intergenicStart","intergenicEnd","overlapLength"))
EMintergenicOverlapsTotal=sum(EMintergenicOverlaps$overlapLength)


###

EMpctExonCov=EMexonOverlapsTotal/exonsTotal
EMpctIntronsCov=EMintronOverlapsTotal/intronsTotal
EMpctintergenicCov=EMintergenicOverlapsTotal/intergenicTotal


########################## E. tales To Melpomene ########################
TMexonOverlaps=read.csv("AlignmentStats/alignmentDepth_EtalToHmelRef_exons.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","exonScaf","exonLabel","exonLabel2","exonStart","exonEnd","dot","dash","dot2","ID","overlapLength"))
TMexonOverlapsTotal=sum(TMexonOverlaps$overlapLength)

####

TMintronOverlaps=read.csv("AlignmentStats/alignmentDepth_EtalToHmelRef_introns.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intronScaf","intronLabel","intronLabel2","intronStart","intronEnd","dot","dash","dot2","ID","overlapLength"))
TMintronOverlapsTotal=sum(TMintronOverlaps$overlapLength)
###

TMintergenicOverlaps=read.csv("AlignmentStats/alignmentDepth_EtalToHmelRef_intergenic.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intergenicScaf","intergenicStart","intergenicEnd","overlapLength"))
TMintergenicOverlapsTotal=sum(TMintergenicOverlaps$overlapLength)


###

TMpctExonCov=TMexonOverlapsTotal/exonsTotal
TMpctIntronsCov=TMintronOverlapsTotal/intronsTotal
TMpctintergenicCov=TMintergenicOverlapsTotal/intergenicTotal


########################## Bombyx To Melpomene ########################
BMexonOverlaps=read.csv("AlignmentStats/alignmentDepth_BmorToHmel_exons.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","exonScaf","exonLabel","exonLabel2","exonStart","exonEnd","dot","dash","dot2","ID","overlapLength"))
BMexonOverlapsTotal=sum(BMexonOverlaps$overlapLength)

####

BMintronOverlaps=read.csv("AlignmentStats/alignmentDepth_BmorToHmel_introns.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intronScaf","intronLabel","intronLabel2","intronStart","intronEnd","dot","dash","dot2","ID","overlapLength"))
BMintronOverlapsTotal=sum(BMintronOverlaps$overlapLength)
###

BMintergenicOverlaps=read.csv("AlignmentStats/alignmentDepth_BmorToHmel_intergenic.bed",sep="\t", header = F, col.names = c("alignedScaf","alignedStart","alignedEnd","intergenicScaf","intergenicStart","intergenicEnd","overlapLength"))
BMintergenicOverlapsTotal=sum(BMintergenicOverlaps$overlapLength)


###

BMpctExonCov=BMexonOverlapsTotal/exonsTotal
BMpctIntronsCov=BMintronOverlapsTotal/intronsTotal
BMpctintergenicCov=BMintergenicOverlapsTotal/intergenicTotal

#####
df <- data.frame(labels=c(rep(c("exons","introns","intergenic"),4)), 
                 align=c(rep("H. melpomene DISCO to H. melpomene ref",3),rep("H. erato ref to H. melpomene ref",3),rep("E. tales to H. melpomene ref",3),rep("B. mori to H. melpomene ref",3)), 
                 values=c(pctExonCov,pctIntronsCov,pctintergenicCov,EMpctExonCov,EMpctIntronsCov,EMpctintergenicCov,TMpctExonCov,TMpctIntronsCov,TMpctintergenicCov,BMpctExonCov,BMpctIntronsCov,BMpctintergenicCov))

df$labels=factor(df$labels, levels=c("exons","introns","intergenic"))
df$align=factor(df$align, levels=c("H. melpomene DISCO to H. melpomene ref","H. erato ref to H. melpomene ref","E. tales to H. melpomene ref","B. mori to H. melpomene ref"))

x <- ggplot(df)+
  geom_bar(aes(x=labels,y=values,fill=align), stat="identity",position = "dodge")+
  geom_text(aes(x=labels,y=values+.01,label=round(values,4))) +
  labs(x="",y="Fraction",title="Fraction Aligned")
x

