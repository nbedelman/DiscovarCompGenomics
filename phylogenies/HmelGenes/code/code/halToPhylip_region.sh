#!/bin/bash

#halToMaf.slurm



refGenome=$1
halFile=$2
inFile=$3
targetGenomes=$4
mafFile=mafs/$(basename $inFile .bed).maf
phylipFile=phylips/$(basename $inFile .bed).fa


hal2maf --noAncestors --noDupes --onlyOrthologs --targetGenomes $targetGenomes --refGenome $refGenome --refTargets $inFile $halFile $mafFile
code/mafToPhylip.py $mafFile $phylipFile
