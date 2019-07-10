#!/bin/bash

#SBATCH -J getTrees
#SBATCH -p serial_requeue,shared,general
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 0-03:00
#SBATCH --mem=1000
#SBATCH -e logs/getTrees.err
#SBATCH -o logs/getTrees.out

genoDir=$1
phymlPath=$2
genoFile=$3
treeDir=$4
windSize=$5

python $genoDir/phylo/phyml_sliding_windows.py -w $windSize -g $genoFile --outgroup 4 \
--phyml $phymlPath --optimise tlr -p $treeDir/$(basename $genoFile .geno) \
--bootstrap 100 --individuals 1,2,3,4 -T 32

##this outputs a gz file, but QuIBL takes just a txt, so unzip it.
gunzip $treeDir/$(basename $genoFile .geno).trees.gz
