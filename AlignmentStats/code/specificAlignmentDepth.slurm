#!/bin/bash

#getAlignmentDepth.slurm

#SBATCH -n 16
#SBATCH -t 0-05:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o getAlignmentDepth.out
#SBATCH -e getAlignmentDepth.err

module load hal

halPath=$1
referenceGenome=$2
outputWig=$3
outputBed=$4
targetGenome=$5
depth=$6

halAlignmentDepth --outWiggle $outputWig --targetGenomes $targetGenome --noAncestors $halPath $referenceGenome
wigToBed.py $outputWig $outputBed $depth 0 1
