#!/bin/bash

genoDir=$1
genoFile=$2
treeDir=$3
phymlPath=$4
quiblTemplate=$5
quiblPath=$6

python $genoDir/phylo/phyml_sliding_windows.py -w 5000 -g $genoFile --outgroup 4 \
--phyml $phymlPath --optimise tlr -p $treeDir/$(basename $genoFile .geno) --bootstrap 100 \
--individuals 1,2,3,4 -T 32

##this outputs a gz file, but QuIBL takes just a txt, so unzip it.
gunzip $treeDir/$(basename $genoFile .geno).trees.gz

code/bootstrapQuIBL.slurm $treeDir/$(basename $genoFile .geno).trees $quiblTemplate 4 \
$quiblPath  $(basename $genoFile .geno).quibl.out 100 500 $(basename $genoFile .geno)
