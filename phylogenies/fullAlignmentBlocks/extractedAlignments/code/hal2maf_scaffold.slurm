#!/bin/bash

#SBATCH -J hal2maf_scaffold
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 1-00:00
#SBATCH -p general
#SBATCH -o hal2maf_scaffold.out
#SBATCH -e hal2maf_scaffold.err


module load hal

halFile=$1
refGenome=$2
scafName=$3
mafFile=$4

hal2maf --noAncestors --refGenome $refGenome --refSequence $scafName $halFile $mafFile
