#!/bin/bash

#SBATCH -J multiTreesPhyML
#SBATCH -p general,shared,serial_requeue
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 2-00:00
#SBATCH --mem=4000
#SBATCH -e multiTrees.err
#SBATCH -o multiTrees.out

#Computing trees in windows using iqtree

module load centos6/0.0.1-fasrc01  gcc/7.1.0-fasrc01 iqtree/1.5.5-fasrc01
module load mafft

fastaDir=$1
treeDir=$2
end=$3
#rootTaxon=$6


for file in $(ls $fastaDir|head -$end|tail -50)
do
  if ! [ -f $treeDir/$(basename $file .fa).treefile ];then
mafft --quiet --auto --inputorder $fastaDir/$file > $fastaDir/$file.realigned.tmp
iqtree-omp -quiet -s $fastaDir/$file.realigned.tmp -m MFP -pre $treeDir/$(basename $file .fa) -safe -nt AUTO #-o $rootTaxon
rm -f $fastaDir/$file.realigned.tmp
fi
done

#python code/findCommonTrees.py $outFile $(basename $outFile .tsv).foundTrees.csv
