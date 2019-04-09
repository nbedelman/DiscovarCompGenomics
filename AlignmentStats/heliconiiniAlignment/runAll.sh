#!/bin/bash

#get alignment coverage so I can find full alignment blocks for phylogenies.

mkdir data
mkdir code

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/heliconiines.hal data
#cp /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/AlignmentStats/getAlignmentDepth.slurm code
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/scripts/wigToBed.py code

export PATH=$PATH:code/

#manually changed the value of how many genomes should align in getAlignmentDepth.slurm from 24 to 23

sbatch code/getAlignmentDepth.slurm data/heliconiines.hal Heliconius_melpomene_melpomene_hmel2 heliconiines_fullAlignmentBlocks_toHmelRef.wig heliconiines_fullAlignmentBlocks_toHmelRef.bed
