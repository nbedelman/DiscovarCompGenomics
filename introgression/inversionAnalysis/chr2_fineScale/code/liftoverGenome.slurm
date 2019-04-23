#!/bin/bash

#liftoverGenome.slurm


#SBATCH -n 1
#SBATCH -t 0-10:00
#SBATCH --mem=5000
#SBATCH -p general
#SBATCH -o liftoverGenome.out
#SBATCH -e liftoverGenome.err

halFile=$1
refGenome=$2
refMap=$3
newRefGenome=$4


#windows=$(awk '{print $4}' $refMap|sort|uniq)
#for window in $windows
#do awk -v win=$window '$4==win' $refMap | \
#halLiftover --outPSL $halFile $refGenome /dev/stdin $newRefGenome $newRefGenome.$refMap.psl |\
halLiftover --outPSL $halFile $refGenome $refMap $newRefGenome $newRefGenome.$refMap.psl |\
python code/parsePSL.py $newRefGenome.$refMap.psl 100 > $newRefGenome.$refMap
#awk -v win=$window '$2>0 {print $4"\t"$7"\t"$8"\t"win"\t.\t"$9}' $newRefGenome.$(basename $refMap .bed).bed >> $newRefGenome.$refMap
done
