#!/bin/bash

seqDir=$1
genoDir=$2
GENODir=$3
windSize=$4
overallBase=$5
combinedFile=$6

for phyFile in $(ls $seqDir/*)
do base=$(basename $phyFile .phy)
python $genoDir/seqToGeno.py -s $phyFile -g $GENODir/$base.geno -f phylip -C $base  -S 1 2 3 4
grep -v "#" $GENODir/$base.geno  >> $combinedFile
done

python $genoDir/ABBABABAwindows.py -w $windSize -g $combinedFile -f haplo -P1 1 1 -P2 2 2 -P3 3 3  -O 4 4 -T 16 \
-o $overallBase.popStats.csv --haploid 1,2,3,4

code/bootstrapFd.sh $overallBase.popStats.csv 100 500 $overallBase.popStats.boots.out

python code/getABBootStats.py $overallBase.popStats.boots.out $overallBase.popStats.boots.summary
