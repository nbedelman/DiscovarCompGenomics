#!/bin/bash

#SBATCH -J extractChromMaf
#SBATCH -p general
#SBATCH -t 0-5:00
#SBATCH -n 1
#SBATCH --mem=5000
#SBATCH -e extractChromMaf.err
#SBATCH -o extractChromMaf.out

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

halFile=$1
refGenome=$2
refSequenceList=$3

for refSequence in $refSequenceList
do hal2maf --refGenome $refGenome --refSequence $refSequence --onlyOrthologs --noAncestors $halFile $refSequence\_to$refGenome.maf
done
