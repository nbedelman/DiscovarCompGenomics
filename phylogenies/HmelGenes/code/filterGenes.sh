#!/bin/bash

#if the gene is aligned among at least 22 species, realign it and then make sure there's less than 10% missing data
fastaFile=$1
outDir=$2

numSpecies=$(grep ">" $fastaFile|wc -l |awk '{print $1}')
noGap=$(msa_view --gap-strip ANY $fastaFile|bioawk -c fastx '{print length($seq)}' |head -1)
totalLength=$(bioawk -c fastx '{print length($seq)}' $fastaFile|head -1)
