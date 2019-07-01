#!/bin/bash

#SBATCH -J runQuIBL
#SBATCH -n 1
#SBATCH -p serial_requeue,shared,general
#SBATCH -t 0-10:00
#SBATCH --mem=5000
#SBATCH -e logs/runQuIBL.err
#SBATCH -o logs/runQuIBL.out

treeFile=$1
quiblTemplate=$2
outgroup=$3
quiblLoc=$4
quiblSummary=$5
numBoots=$6
numLoci=$7
baseID=$8

echo "C,lambda,introProp,dBIC" > $quiblSummary
for boot in $(seq 1 $numBoots)
do
quiblInput=quibl/$baseID.quiblInput.$boot.txt
quiblOut=quibl/$baseID.quiblOut.$boot.csv
treeSub=quibl/$baseID.quiblBoot.$boot.nwk
shuf -n $numLoci -r $treeFile > $treeSub
sed "s:TREEFILE:$treeSub:" $quiblTemplate | sed "s:OUTGROUP:$outgroup:" | sed "s:OUTPATH:$quiblOut:" > $quiblInput
code/runQuIBL.slurm $quiblLoc $quiblInput $quiblOut $quiblSummary $numLoci 
done
