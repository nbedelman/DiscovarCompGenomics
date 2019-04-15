#!/bin/sh

#mafDirToFasta.slurm

#SBATCH -n 1
#SBATCH -t 0-02:00
#SBATCH --mem=1000
#SBATCH -p general
#SBATCH -o mafDirToFasta.out
#SBATCH -e mafDirToFasta.err


mafDir=$1

for maf in "$mafDir/*"
do fasta=$(basename $maf .maf).fa
msa_view $maf > $fasta
done
