#!/bin/bash

#threadTrees.slurm

#SBATCH -J threadTrees
#SBATCH -p general
#SBATCH -n 1
#SBATCH -t 0-05:00
#SBATCH --mem=1000
#SBATCH -e threadTrees.err
#SBATCH -o threadTrees.out

chromNumber=$1


for f in $(ls fastas|grep chr$chromNumber\_)
do cp fastas/$f chr$chromNumber\_fastas
done

cd chr$chromNumber\_fastas
sbatch ../code/multiTrees.slurm /n/regal/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/introgression/eratoClade_10KBwindows/trees allTrees_$chromNumber.tsv
cd ..
