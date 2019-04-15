#!/bin/bash

#slurm wrapper for faDirToNexus.py

#SBATCH -J faDirToNexus
#SBATCH -p general
#SBATCH -n 1
#SBATCH --mem=1000
#SBATCH -t 00-01:00
#SBATCH -e faDirToNexus.err
#SBATCH -o faDirToNexus.out

faDir=$1
output=$2
numIntros=$3

code/faDirToNexus.py $faDir $output $numIntros
