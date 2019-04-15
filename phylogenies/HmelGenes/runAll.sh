#!/bin/bash

#SBATCH -J h2m_codingBlocks                  # A single job name for the array
#SBATCH -n 1                       # Number of cores
#SBATCH -p general,shared,serial_requeue         # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 1-00:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o h2m_codingBlocks.out
#SBATCH -e h2m_codingBlocks.err

mkdir -p data
mkdir -p beds
mkdir -p mafs
mkdir -p fastas
mkdir -p ORFs

######## Setup environment ###########

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

module load ucsc

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/heliconius_melpomene_melpomene_hmel2_core_32_85_1.gff  data/HmelRef.gff
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal


######### Define Variables
origGffFile=data/HmelRef.gff
gffFile=data/HmelRef.header.gff
CDSgffFile=data/HmelRef.header.CDSAsExon.gff
genePred=data/HmelRef.header.CDSAsExon.gp
bedFile=data/fullGenes.header.CDSAsExon.bed

halFile=data/fullAlignment.hal
simpBedFile=data/fullGenes.oneTrans.bed
outBase=fullGeneModel

refGenome=HmelRef

end=$1

########## Run Code #########
#First, get gene pred of ONLY coding sequence. This means changing the IDs of CDS to exon and getting rid of exons.
echo '##gff-version 3' > $gffFile
cat $origGffFile >> $gffFile
awk '$3!="exon"' $gffFile |sed 's/CDS/exon/g' > $CDSgffFile

##Then, convert gff > genePred > bed.
gff3ToGenePred $CDSgffFile $genePred
genePredToBed $genePred  $bedFile

just take the first transcript for each gene
grep t1  $bedFile > $simpBedFile

code/batchGenes.sh $simpBedFile $outBase

for gene in $(ls beds|head -$end |tail -1000)
do if ! [ -f fastas/$(basename $gene .bed).fa ];then
code/halToMaf_codingRegion.sh $refGenome $halFile beds/$gene
fi
done
