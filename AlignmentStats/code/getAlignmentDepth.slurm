#!/bin/bash

#getAlignmentDepth.slurm

#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o getAlignmentDepth.out
#SBATCH -e getAlignmentDepth.err

module load hal

halPath=$1
referenceGenome=$2
outputWig=$3
outputBed=$4

halAlignmentDepth --outWiggle $outputWig --noAncestors $halPath $referenceGenome
wigToBed.py $outputWig $outputBed 24 0 1
