#!/bin/bash

#SBATCH -J convertAndMergePhylips
#SBATCH -p serial_requeue,shared,general
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 0-08:00
#SBATCH --mem=1000
#SBATCH -e logs/convertAndMerge.err
#SBATCH -o logs/convertAndMerge.out

genoDir=$1
phyDir=$2
GENODir=$3
combinedFile=$4

header=$(echo \#CHROM POS 1 2 3 4)
echo $header > $combinedFile

for phyFile in $(ls $phyDir/*)
do base=$(basename $phyFile .dat)
python $genoDir/seqToGeno.py -s $phyFile -g $GENODir/$base.geno -f phylip -C $base  -S 1 2 3 4
grep -v "#" $GENODir/$base.geno  >> $combinedFile
done
