#!/bin/bash

### This is a runAll script to construct windows of a chosen size in a chosen species that correspond to the H. erato demaphoon reference order. ###
## 10KB windows every 50KB
############################################
################ Setup Environment #########
############################################

mkdir data
mkdir code

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal
ln -s /n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/Herato2.fasta data/refGenome.fa

module load bedtools2

chmod u+x code/*
export PATH=code/:$PATH

############################################
################ Define Variables ##########
############################################

#make tsv
refFasta=data/refGenome.fa
transitionTSV=data/Herato.transitions.tsv

#make windows

refBed=HeratoWindows.bed
refGenome=HeraRef
windowSize=10000
slideSize=10000



############################################
################ Run Code ##################
############################################

#Make erato transition tsv file
code/makeEratoTransitionTSV.py $refFasta $transitionTSV

#Get the HmelRef map directly from the transition tsv file
awk '{print $7"\t"$8"\t"$9}' $transitionTSV > $refGenome.map.bed

#calculate appropriate windows
sbatch code/makeReferenceWindows.slurm $transitionTSV $refBed $windowSize $slideSize
