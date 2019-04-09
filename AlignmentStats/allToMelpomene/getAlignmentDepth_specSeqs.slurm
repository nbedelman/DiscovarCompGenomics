#!/bin/bash

#getAlignmentDepth_specSeqs.slurm

#SBATCH -n 16
#SBATCH -t 0-10:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o getAlignmentDepth_specSeqs.out
#SBATCH -e getAlignmentDepth_specSeqs.err

halPath=$1
referenceGenome=$2
outputWig=$3
outputDir=$4
seq=$5

mkdir -p $outputDir

halAlignmentDepth --outWiggle $outputDir/$outputWig\_$seq --refSequence $seq --noAncestors $halPath $referenceGenome
