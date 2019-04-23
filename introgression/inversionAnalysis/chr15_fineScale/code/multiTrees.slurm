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

for file in chr*.fa
do raxmlHPC -s $file -n $file -m GTRCAT -p 12345 --silent -w $outDir
echo $(basename $file .fa) $(cat ../trees/RAxML_bestTree.$file) >> ../allTrees.tsv
done

cd ..
python code/filterTrees.py allTrees.tsv filteredTrees.csv
