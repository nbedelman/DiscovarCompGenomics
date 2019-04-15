#!/bin/bash

#getABBABABA.slurm
#for melmel-cydtim-parnum

#SBATCH -J ABBABABA                # A single job name for the array
#SBATCH -n 1                       # Number of cores
#SBATCH -p general         # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-05:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o ABBABABA.out        # Standard output
#SBATCH -e ABBABABA.err        # Standard error

fastaDir=$1
summaryCSV=$2

source activate my_root

for fasta in $fastaDir/*
do
sed 's/ //g' $fasta > $fasta\_tmp
vcfOut=VCFs/$(basename $fasta\_tmp .fa_tmp).vcf
snp-sites -v -o $vcfOut $fasta\_tmp
alnLength=$(head -2 $vcfOut|tail -1|cut -d "=" -f 4|cut -d ">" -f 1)
numSNPs=$(grep -v "#" $vcfOut|wc -l)
numBial=$(grep -v "#" $vcfOut|grep -v "," |wc -l)
numABBA=$(grep -v "#" $vcfOut|grep -v "," |awk '$10==$13'|awk '$10!=$11'|awk '$11==$12'|wc -l)
numBABA=$(grep -v "#" $vcfOut|grep -v "," |awk '$10==$12'|awk '$10!=$11'|awk '$11==$13'|wc -l)
echo $(basename $fasta .fa),$alnLength,$numSNPs,$numBial,$numABBA,$numBABA >> $summaryCSV
rm -f $fasta\_tmp
done
