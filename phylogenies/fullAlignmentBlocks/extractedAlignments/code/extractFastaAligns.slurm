#!/bin/sh

#SBATCH -J extractFastaAligns
#SBATCH -p general,shared
#SBATCH -n 1
#SBATCH -t 0-10:00
#SBATCH --mem=4000
#SBATCH -e extractFastaAligns.err
#SBATCH -o extractFastaAligns.out

directory=$1
mafDir=$2
fastaDir=$3
targetGenomes=$4
refGenome=$5
halFile=$6

for gene in $(ls "$directory/"*)
do code/halToMaf_codingRegion.sh $refGenome $halFile $gene $mafDir $fastaDir $targetGenomes
done
