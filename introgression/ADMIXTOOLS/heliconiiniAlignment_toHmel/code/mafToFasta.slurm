#!/bin/bash

#SBATCH -J mafToFasta
#SBATCH -p general
#SBATCH -t 0-10:00
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -e mafToFasta.err
#SBATCH -o mafToFasta.out

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

mafDir=$1

# for maf in $mafDir/*
# do fasta=fastas/$(basename $maf .maf).fa
# msa_view --in-format MAF --out-format FASTA -G 1 -m $maf > $fasta
# done

emptyFastas=$(find  fastas/ -type f -empty)
for fasta in $emptyFastas
do maf=mafs/$(basename $fasta .fa).maf
code/singleMafToFasta.slurm $maf
done
