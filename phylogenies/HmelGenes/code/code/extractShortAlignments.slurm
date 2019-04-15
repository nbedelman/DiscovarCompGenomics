#!/bin/sh

#SBATCH -J extractShortAlignments
#SBATCH -p serial_requeue,shared,general
#SBATCH -t 1-00:00
#SBATCH --mem=5000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -e extractShortAlignments.err
#SBATCH -o extractShortAlignments.out

scafMapFolder=$1
searchSuffix=$2

for mafFile in $(ls $scafMapFolder/|grep $searchSuffix)
do mafSpeciesSubset $scafMapFolder/$mafFile targetSeqs.txt tmp1.maf
mafFilter -minRow=7 -maxRow=7 -minCol=200 tmp1.maf > tmp2.maf
code/getSingleCopy.py tmp2.maf
code/getFullAlignMaf.py tmp2_singleCopy.maf .1 tmp3_singleCopy
mafSplit tmp3_singleCopy_fullAlign.scafMap.bed mafs/ tmp3_singleCopy_fullAlign.maf
rm -f tmp*
done


for i in mafs/*
do
msa_view $i > fastas/$(basename $i .maf).fa
done
