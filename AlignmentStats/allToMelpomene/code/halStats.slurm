#!/bin/bash

#this code runs halStats with the percentID flag so we can see the percent identity between Hmel and all other genomes.
#halStats.slurm

#SBATCH -n 16
#SBATCH -t 0-10:00
#SBATCH --mem=50000
#SBATCH -p general
#SBATCH -o halStats.out
#SBATCH -e halStats.err

halFile=$1
refGenome=$2

halStats --percentID $refGenome $halFile > $(basename $halFile .hal)\_percentID.stats
