#!/bin/bash

#SBATCH -J getFourPopStats
#SBATCH -p serial_requeue,shared,general
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 0-01:00
#SBATCH --mem=1000
#SBATCH -e logs/getFourPopStats.err
#SBATCH -o logs/getFourPopStats.out

genoDir=$1
genoFile=$2
windSize=$3
outFile=$4
summaryFile=$5
summaryBoot=$6
numBoots=$7
numLoci=$8

python $genoDir/ABBABABAwindows.py -w $windSize -g $genoFile -f haplo -P2 2 2 -P1 1 1 -P3 3 3 -O 4 4 -T 16 --haploid 1,2,3,4 -o $outFile
python code/getSummaryStats.py $outFile $summaryFile
code/bootstrapFd.sh $outFile $numBoots $numLoci $outFile.boots
python code/getABBootStats.py $outFile.boots $summaryBoot
