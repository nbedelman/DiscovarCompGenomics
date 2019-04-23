#!/bin/bash

#SBATCH -J subsetAndSingleCopy
#SBATCH -n 1
#SBATCH --mem=5000
#SBATCH -t 0-02:00
#SBATCH -p general
#SBATCH -o subsetAndSingleCopy.out
#SBATCH -e subsetAndSingleCopy.err

exclusions=$1
mafFile=$2
list=$3

mafBase=$(echo $mafFile|cut -d "." -f 1 )
subsetMaf=$mafBase.$list.maf

grep -Ev "$exclusions" $mafFile > $subsetMaf

code/getSingleCopy.py $subsetMaf
