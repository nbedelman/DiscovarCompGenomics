#!/bin/bash

#extractMaf.sh
#Nate Edelman, 6/29/16
#Takes a hal file and a reference genome
#returns a maf file of the alignment

#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH --mem=50000
#SBATCH -p serial_requeue
#SBATCH -o halToMaf.out
#SBATCH -e halToMaf.err

halFile=$1
refGenome=$2
mafFile=$(basename $halFile .hal)_/$scafName/.maf

hal2maf --noAncestors --refGenome $refGenome $halFile $mafFile
