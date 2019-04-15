#!/bin/bash

#SBATCH -J multiTreesPhyML
#SBATCH -p general
#SBATCH -n 1
#SBATCH -t 2-00:00
#SBATCH --mem=4000
#SBATCH -e multiTrees.err
#SBATCH -o multiTrees.out

#Computing trees in windows for michael's analysis.

module load phyml/2014Oct16-fasrc01

fastaDir=$1
phylipDir=$2
outFile=$3
minSize=$4
maxSize=$5

for file in $(ls $fastaDir)
do len=$(bioawk -c fastx '{print length($seq)}' $fastaDir/$file|head -1)
if [[ $len -ge $minSize && $len -le $maxSize ]];then
phyFile=$phylipDir/$(basename $file .fa).phylip
msa_view $fastaDir/$file --out-format PHYLIP > $phyFile 
phyml -i $phyFile -d nt -n 1 -s BEST -o 'tlr' -b 10000 --n_rand_starts 10 --quiet -q
echo $(basename $file .fa) $(cat $phyFile\_phyml_tree) $(grep "Log-likelihood" $phyFile\_phyml_stats|awk '{print $3}')>> $outFile
fi
done

python code/rootTrees.py $outFile $(basename $outFile .tsv).rooted.tsv
awk '{print $2}' $(basename $outFile .tsv).rooted.tsv |grep -v NA > $(basename $outFile .tsv).rooted.newickOnly.txt

python code/findCommonTrees.py $outFile $(basename $outFile .tsv).foundTrees.csv
