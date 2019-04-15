#!/bin/bash

#halToMaf.slurm

#SBATCH -J h2m_codingBlocks                  # A single job name for the array
#SBATCH -n 16                       # Number of cores
#SBATCH -p serial_requeue          # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 0-0:45                  # Maximum execution time (D-HH:MM)
#SBATCH -o h2m_codingBlocks_%A_%a.out        # Standard output
#SBATCH -e h2m_codingBlocks_%A_%a.err        # Standard error

refGenome=$1
halFile=$2
outBase=$3
inFile=beds/$outBase\_"${SLURM_ARRAY_TASK_ID}"_HMEL*.bed
mafFile=mafs/$(basename $inFile .bed).maf
fastaFile=fastas/$(basename $inFile .bed).fa
strand=$(awk '{print $6}' $inFile)


hal2maf --noAncestors --noDupes --refGenome $refGenome --refTargets $inFile $halFile $mafFile
code/mafToFasta.py $mafFile $fastaFile
if [ $strand == '-' ];then
msa_view -V $fastaFile > $fastaFile.tmp && mv -f $fastaFile.tmp $fastaFile
fi
