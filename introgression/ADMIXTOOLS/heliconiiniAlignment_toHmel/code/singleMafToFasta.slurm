#!/bin/bash

#SBATCH -J mafToFasta
#SBATCH -p general
#SBATCH -t 0-1:00
#SBATCH -n 1
#SBATCH --mem=5000
#SBATCH -e singleMafToFasta.err
#SBATCH -o singleMafToFasta.out

export PATH=$PATH:/n/mallet_lab/edelman/software/kentUtils/bin
module load libpng/1.5.21-fasrc01

mafFile=$1
scaf=$(basename $mafFile _toHeliconius_melpomene_melpomene_hmel2.maf)

length=$(awk '{print $6}' $mafFile |head -5|tail -1)
start=0
end=100001
while [ $end -lt $length ]
do echo $scaf $start $end >> $scaf.bed
let start="$end+1"
let end="$start+100000"
done
echo $scaf $start $length >> $scaf.bed

mafSplit $scaf.bed . $mafFile

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

for i in $scaf.*.maf
do fasta=fastas/$(basename $i .maf).fa
msa_view --in-format MAF --out-format FASTA -G 1 -m $i > $fasta
done

mv $scaf.*.maf mafs
