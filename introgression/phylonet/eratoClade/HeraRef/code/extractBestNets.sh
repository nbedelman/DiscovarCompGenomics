#!/bin/bash

inFile=$1
outBase=$2
outNoWeight=$outBase.noWeight.tree
outOnlyBest=$outBase.onlyBest.tree

grep Rank $inFile|cut -d ";" -f 4|cut -d ":" -f 2- |sed 's/$/;/g'>> $outBase.tree
grep Rank $inFile|cut -d ";" -f 4|cut -d ":" -f 2- | sed 's/::[0-1].[0-9]*//g'|sed 's/$/;/g'>>$outNoWeight
grep "Rank = 0" $inFile|cut -d ";" -f 4|cut -d ":" -f 2- | sed 's/::[0-1].[0-9]*//g'|sed 's/$/;/g'>>$outOnlyBest
