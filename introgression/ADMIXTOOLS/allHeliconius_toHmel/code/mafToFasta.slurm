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
fastaDir=$2

for maf in $mafDir/*
do fasta=$fastaDir/$(basename $maf .maf).fa
msa_view --in-format MAF --out-format FASTA -G 1 -m $maf > $fasta
done
