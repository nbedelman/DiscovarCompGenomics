#!/bin/bash

#slurm wrapper for vcfToDistances.py

#SBATCH -J vcfToDistances
#SBATCH -p general
#SBATCH -n 1
#SBATCH -t 0-10:00
#SBATCH --mem=4000
#SBATCH -e vcfToDistances.err
#SBATCH -o vcfToDistances.out

VCFDir=$1
outFile=$2

python code/vcfToDistances.py $VCFDir $outFile
