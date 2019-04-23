#!/bin/bash

#Making this runAll into a slurm file because each script takes too little time. May regret this decision, but can always revert back to sending small slurm jobs.

#SBATCH -J getSingleCopyEverything
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 3-00:00
#SBATCH -p general
#SBATCH -o getSingleCopyEverything.out
#SBATCH -e getSingleCopyEverything.err

#script to get the "usable parts" of the hal alignment. This means projecting the alignment onto good genomes and only taking the single-copy regions.

######### Setup Environment ##############
mkdir -p data
mkdir -p code
mkdir -p mafs

source /n/sw/progressiveCactus-latest/progressiveCactus/environment
export PATH=$PATH:/n/mallet_lab/edelman/software/kentUtils/bin
module load libpng/1.5.21-fasrc01
module load legacy/0.0.1-fasrc01 centos6/openssl-1.0.1f
module load bedtools


ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal
ln -s /n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Herato2.fasta_1kb.fa data/HeraRef.fa
ln -s /n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Hmel2.fa_1kb.fa data/HmelRef.fa
ln -s /n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/DAS_09-132_e_tales_a-scaffolds.fasta_1kb.fa data/Etal.fa
ln -s /n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Bombyx_mori_ASM15162v1_-_scaffolds.fa_1kb.fa data/Bmor.fa
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/heliconius_melpomene_melpomene_hmel2_core_32_85_1.gff data/Hmel.gff


########## Define Variables #############
refGenomes=HmelRef
refGenomes=$(echo $refGenomes|sed 's/,/\n/g')
fullRefGenomes=Bmor,HmelRef,HeraRef
fullRefGenomes=$(echo $fullRefGenomes|sed 's/,/\n/g')
helRefGenomes=Etal,HmelRef,HeraRef
helRefGenomes=$(echo $helRefGenomes|sed 's/,/\n/g')

halFile=data/fullAlignment.hal
gffFile=data/Hmel.gff
CDSBed=data/HmelRef.CDS.bed

#the exclusion lists will be for things like introgression analysis. I have several subsets here to use depending on the analysis desired.
allHel="Pxyl|Bmor|Lacc|Ppol|Dple|Bany|Mcin"
#helSlim doesn't have the hybrid, and only includes one of the melpomene and erato genomes.
helSlim="Pxyl|Bmor|Lacc|Ppol|Dple|Bany|Mcin|HeraDisco|HmelDisco|HeraHhimHyb"
melClade="Pxyl|Bmor|Lacc|Ppol|Dple|Bany|Mcin|Avan|HeraDisco|HmelDisco|HeraHhimHyb|Hdem|Hsar|Htel|Hhsa|Hhim|HeraRef"
eraClade="Pxyl|Bmor|Lacc|Ppol|Dple|Bany|Mcin|Avan|Ldor|HmelRef|Hcyd|Htim|Hbes|Hnum|Hpar|HeraDisco|HmelDisco|HeraHhimHyb"

exclusionList=$allHel,$helSlim,$melClade,$eraClade
exclusionList=$(echo $exclusionList|sed 's/,/\n/g')
######### run code ##########
#extract full alignment for each of the refGenomes. Do this by scaffold for efficiency.
#then, use grep to get sub-alignments for each species set
#and extract only the single-copy alignment blocks
for r in $refGenomes
do
  mkdir -p $r\_mafs
  scaffolds=$(faSize -detailed data/$r.fa |awk '{print $1}')
  for scaffold in $scaffolds
  do
    mafFile=$r\_mafs/$r\_$scaffold.full.maf
    # extractFullMaf=`sbatch code/hal2maf_scaffold.slurm $halFile $r $scaffold $mafFile|cut -d " " -f 4`
    # allSingleCopy=`sbatch --dependency=afterok:$extractFullMaf code/subsetAndSingleCopy.slurm "None" $mafFile "full"`
    code/hal2maf_scaffold.slurm $halFile $r $scaffold $mafFile
    code/subsetAndSingleCopy.slurm "None" $mafFile "all"
    for ex in $exclusionList
    do #figure out which list we have
      if [ $ex == $allHel ];then
        list=allHel
      elif [ $ex == $helSlim ];then
        list=helSlim
      elif [ $ex == $melClade ];then
        list=melClade
      else
        list=eraClade
      fi
      # singleCopy=`sbatch --dependency=afterok:$extractFullMaf code/subsetAndSingleCopy.slurm $ex $mafFile $list`
      code/subsetAndSingleCopy.slurm $ex $mafFile $list
    done
  done
done

#after I ran this code, I realized that I needed to re-do the fully-aligned sites analysis because I originally did not consider whether regions were duplicated or not.
mkdir $r\_beds
mkdir $r\_fastas

grep -v "#" $gffFile| awk '$3=="CDS" {print $1"\t"$4"\t"$5"\t"$9"\t.\t"$7}' > $CDSBed

for i in $r\_mafs/*.all_singleCopy.maf; do code/getFullyAligned.py $i 25 $r\_mafs/$(basename $i .maf); done
cat $r\_mafs/*.fullAlign*bed > $r\_mafs/HmelRef.fullyAlignedSites.bed
awk '$3-$2 > 150' $r\_mafs/HmelRef.fullyAlignedSites.bed  > $r\_mafs/HmelRef.fullyAlignedSites.large.bed

bedtools intersect -a $r\_mafs/HmelRef.fullyAlignedSites.large.bed -b $CDSBed |sort -u > $r\_beds/HmelRef.fullyAlignedSites.large.coding.bed
bedtools intersect -v -a $r\_mafs/HmelRef.fullyAlignedSites.large.bed -b $CDSBed > $r\_beds/HmelRef.fullyAlignedSites.large.noncoding.bed

#put each region into its own bed file

code/batchGenes.sh $r\_beds/HmelRef.fullyAlignedSites.large.coding.bed $r\_beds/fullyAligned.Hmel.coding
code/batchGenes.sh $r\_beds/HmelRef.fullyAlignedSites.large.noncoding.bed $r\_beds/fullyAligned.Hmel.noncoding

for gene in $(ls beds)
do code/halToMaf_codingRegion.sh $refGenome $halFile beds/$gene
done
############## all heliconius
mkdir $r\_allHeliconius_fullyAligned

workDir=$r\_allHeliconius_fullyAligned

for i in $r\_mafs/*.allHel_singleCopy.maf; do code/getFullyAligned.py $i 18 $workDir/$(basename $i .maf); done
cat $workDir/*.fullAlign*bed > $workDir/HmelRef.allHel.fullyAlignedSites.bed
awk '$3-$2 > 150' $workDir/HmelRef.allHel.fullyAlignedSites.bed  > $workDir/HmelRef.allHel.fullyAlignedSites.large.bed

bedtools intersect -a $workDir/HmelRef.allHel.fullyAlignedSites.large.bed -b $CDSBed |sort -u > $workDir/HmelRef.allHel.fullyAlignedSites.coding.bed
bedtools intersect -v -a $workDir/HmelRef.allHel.fullyAlignedSites.large.bed  -b $CDSBed > $workDir/HmelRef.allHel.fullyAlignedSites.noncoding.bed

#put each region into its own bed file

code/batchGenes.sh $r\_beds/HmelRef.fullyAlignedSites.large.coding.bed $r\_beds/fullyAligned.Hmel.coding
code/batchGenes.sh $r\_beds/HmelRef.fullyAlignedSites.large.noncoding.bed $r\_beds/fullyAligned.Hmel.noncoding

for gene in $(ls beds)
do code/halToMaf_codingRegion.sh $refGenome $halFile beds/$gene
done
