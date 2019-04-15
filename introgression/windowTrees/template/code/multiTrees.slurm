#!/bin/bash

#SBATCH -J multiTrees
#SBATCH -p general
#SBATCH -n 1
#SBATCH -t 0-10:00
#SBATCH --mem=4000
#SBATCH -e multiTrees.err
#SBATCH -o multiTrees.out

source activate my_root

outDir=$1
outFile=$2

for file in chr*.fa
do raxmlHPC -s $file -n $file -m GTRCAT -p 12345 --silent -w $outDir
echo $(basename $file .fa) $(cat ../trees/RAxML_bestTree.$file) >> $outFile
done

cd ..
python code/rootTrees.py allTrees.tsv rootedTrees.tsv
awk '{print $2}' rootedTrees.tsv |grep -v NA > rootedTrees.newickOnly.txt

python code/findCommonTrees.py allTrees.tsv foundTrees.csv
