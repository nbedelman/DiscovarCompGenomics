#!/bin/bash

#slurm wrapper for faDirToNexus.random.py

#SBATCH -J faDirToNexus
#SBATCH -p general
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 0-05:00
#SBATCH -e faDirToNexus.err
#SBATCH -o faDirToNexus.out

faDir=$1
output=$2
numIntros=$3
totalAvailable=$4
totalDesired=$5
minSize=$6
additionalOptions=$6

code/faDirToNexus.random.py $faDir $output $numIntros $totalAvailable $totalDesired $minSize "$additionalOptions"
