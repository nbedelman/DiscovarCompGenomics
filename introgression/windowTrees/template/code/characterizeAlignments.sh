#!/bin/bash

#create tsv with region name, length of reference sequence, length of ungapped sequence

alignDir=$1
outFile=$2


for f in $alignDir/*
do
region=$(basename $f .fa)
refLength=$(bioawk -c fastx '{print length($seq)}' $f|head -1)
noGap=$(msa_view --gap-strip ANY $f|bioawk -c fastx '{print length($seq)}' |head -1)
echo $region $refLength $noGap >> $outFile
done
