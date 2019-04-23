#!/bin/bash

#halToMaf.slurm

#SBATCH -J h2m_codingBlocks                  # A single job name for the array
#SBATCH -n 1                       # Number of cores
#SBATCH -p general,serial_requeue,shared         # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 1-00:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o h2m_codingBlocks.out        # Standard output
#SBATCH -e h2m_codingBlocks.err        # Standard error

refGenome=$1
halFile=$2
targetGenomes=$3
refTargets=$4

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

windowNames=$(awk '{print $4}' $refTargets|sort|uniq)


#hal2maf wants the start coord of beds to be base 0
awk '{print $1"\t"$2-1"\t"$3"\t"$4}' $refTargets > $refTargets\_adj_tmp.bed

for window in $windowNames
do awk -v window=$window '$4==window {print}' $refTargets\_adj_tmp.bed > beds/$window.bed
#extract the alignment in maf format from the hal file
#msa_view can only deal with a single scaffold, so have to do them all separately
mafBase=mafs/$window
chroms=$(awk '{print $1}' beds/$window.bed|sort|uniq)
for chrom in $chroms
do grep $chrom beds/$window.bed > beds/$window\_$chrom.bed
hal2maf --refGenome $refGenome --targetGenomes $targetGenomes --refTargets beds/$window\_$chrom.bed $halFile $mafBase\_$chrom.maf

code/getSingleCopy.py $mafBase\_$chrom.maf
maf-sort $mafBase\_$chrom\_singleCopy.maf | tail -n +2 > $mafBase\_$chrom\_singleCopy.sorted.maf

fastaBase=fastas/$window
msa_view $mafBase\_$chrom\_singleCopy.sorted.maf -m --gap-strip 1 > $fastaBase\_$chrom\_singleCopy.fa
if [ ! -s $fastaBase\_$chrom\_singleCopy.fa ]; then rm -f $fastaBase\_$chrom\_singleCopy.fa; fi

done
msa_view --aggregate $targetGenomes  $fastaBase\_*> $fastaBase.fa
rm -f $fastaBase\_*
rm -f beds/$window\_*
rm -f mafs/$window\_*
done

rm  -f $refTargets\_adj_tmp.bed
