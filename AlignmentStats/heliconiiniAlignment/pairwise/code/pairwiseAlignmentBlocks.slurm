#!/bin/bash

#getAlignmentDepth.slurm

#SBATCH -n 16
#SBATCH -t 0-05:00
#SBATCH --mem=10000
#SBATCH -p general
#SBATCH -o getAlignmentDepth.out
#SBATCH -e getAlignmentDepth.err


halPath=$1
referenceGenome=$2
outputWig=$3
outputBed=$4
targetGenome=$5
codingBed=$6
intronBed=$7
geneBed=$8
overallOut=$9

halAlignmentDepth --outWiggle $outputWig --targetGenomes $targetGenome --noAncestors $halPath $referenceGenome
wigToBed.py $outputWig $outputBed 1 0 1
codingBases=$(bedtools intersect -a $codingBed -b $outputBed -wo | awk '{print $10}'| paste -sd+ - | bc)
intronBases=$(bedtools intersect -a $intronBed -b $outputBed -wo | awk '{print $10}'| paste -sd+ - | bc)
intergenicBases=$(bedtools intersect -v -a $f -b data/Hmelgene.bed |awk '{print $3-$2+1}'|paste -sd+ - | bc)

echo $targetGenome $codingBases $intronBases $intergenicBases >> $overallOut
