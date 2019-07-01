#!/bin/bash

overallStatsFile=$1
numBoots=$2
numLoci=$3
outFile=$4


for boot in $(seq 1 $numBoots)
do  shuf -n $numLoci -r $overallStatsFile | \
grep -v nan |awk 'BEGIN { FS = "," } ; $9>=0  {fd += $10; c++}; {D += $9; n++} END {print fd/c","D/n}' >> $outFile
done
