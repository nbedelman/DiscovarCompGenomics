#!/bin/bash

#visualizeResults.slurm
#

#SBATCH -J visualize                # A single job name for the array
#SBATCH -n 1                       # Number of cores
#SBATCH -p general         # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-00:10                  # Maximum execution time (D-HH:MM)
#SBATCH -o visualizeResults.out        # Standard output
#SBATCH -e visualizeResults.err        # Standard error

module load gcc/7.1.0-fasrc01

ABBABABBAfile=$1
treeFile=$2
outBase=$3
mapFile=$4
fullChroms=$5

R CMD BATCH --no-save --no-restore "--args $ABBABABBAfile $treeFile $outBase $mapFile $fullChroms" code/windowAnalysis.R
