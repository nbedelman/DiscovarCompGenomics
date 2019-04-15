#!/bin/bash

#halToMaf.slurm

#SBATCH -J h2m_codingBlocks                  # A single job name for the array
#SBATCH -n 1                       # Number of cores
#SBATCH -p general,shared, serial_requeue          # Partition
#SBATCH --mem 4000                 # Memory request (4Gb)
#SBATCH -t 1-00:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o h2m_codingBlocks.out
#SBATCH -e h2m_codingBlocks.err


refGenome=$1
halFile=$2
inFile=$3
mafFile=mafs/$(basename $inFile .bed).maf
fastaFile=fastas/$(basename $inFile .bed).fa
orfOnly=ORFs/$(basename $inFile .bed).ORFs.fa
strand=$(awk '{print $6}' $inFile)


hal2maf --noAncestors --noDupes --onlyOrthologs --refGenome $refGenome --refTargets $inFile $halFile $mafFile
code/mafToFasta.py $mafFile $fastaFile
if [ $strand == '-' ];then
msa_view -V $fastaFile > $fastaFile.tmp && mv -f $fastaFile.tmp $fastaFile
fi
code/getOpenReadingFrames.py $fastaFile $orfOnly $refGenome
