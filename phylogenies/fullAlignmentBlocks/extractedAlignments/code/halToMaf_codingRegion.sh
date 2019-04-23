#!/bin/bash

#halToMaf.sh

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

refGenome=$1
halFile=$2
inFile=$3
mafDir=$4
fastaDir=$5
targetGenomes=$6

mkdir -p $mafDir
mkdir -p $fastaDir

mafFile=$mafDir/$(basename $inFile .bed).maf
fastaFile=$fastaDir/$(basename $inFile .bed).fa




hal2maf --noAncestors --noDupes --refGenome $refGenome --refTargets $inFile --targetGenomes $targetGenomes $halFile $mafFile
code/mafToFasta.py $mafFile $fastaFile
