#!/bin/sh

#halToMaf.slurm

#SBATCH -n 1                       # Number of cores
#SBATCH -p serial_requeue          # Partition
#SBATCH --mem 5000                 # Memory request (5Gb)
#SBATCH -t 0-2:00                  # Maximum execution time (D-HH:MM)
#SBATCH -o h2m_codingBlocks_%A.out        # Standard output
#SBATCH -e h2m_codingBlocks_%A.err        # Standard error



refGenome=$1
halFile=$2
outBase=$3
type=$4
start=$5
end=$6


for number in $(seq $start $end)
do inFile=$type\Beds/$outBase\_$number.bed
mafFile=$type\Mafs/$(basename $inFile .bed).maf
fastaFile=$type\Fastas/$(basename $inFile .bed).fa
hal2maf --noAncestors --refGenome $refGenome --refTargets $inFile $halFile $mafFile
msa_view $mafFile --in-format MAF --seqs Pxyl,Bmor,Lacc,Ppol,Dple,Bany,Mcin,Avan,Etal,Hdem,Hsar,Htel,Hhsa,HeraHhimHyb,Hhim,HeraDisco,HeraRef,Ldor,HmelDisco,HmelRef,Hcyd,Htim,Hbes,Hnum,Hpar> $fastaFile
done
