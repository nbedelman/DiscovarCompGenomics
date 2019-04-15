#!/bin/bash

#slurm wrapper for faDirToNexus.random.py

#SBATCH -J faDirToNexus
#SBATCH -p general
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 00-05:00
#SBATCH -e faDirToNexus.err
#SBATCH -o faDirToNexus.out

faDir=$1
outBase=$2
numIntros=$3
numFiles=$4
numEntries=$5
minSize=$6
maxSize=$7
additionalOptions=$8

code/faDirToNexus.random.new.py $faDir $outBase $numIntros $numFiles $numEntries $minSize $maxSize "$additionalOptions"
