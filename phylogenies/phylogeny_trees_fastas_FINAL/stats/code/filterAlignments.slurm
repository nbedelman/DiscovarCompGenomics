#!/bin/bash

#SBATCH -J filterAlignments
#SBATCH -p serial_requeue,shared,general
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -e filterAlignments.err
#SBATCH -o filterAlignments.out

faDir=$1
statsFile=$2
treeDir=$3
outDir=$4


#Calculate trees and output relevant stats
for fa in $(ls $faDir)
do
id=$(basename $fa .fas)
outNex=$outDir/$id.fullAlign.nex

# Provide feedback on screen.
#echo -n "Analyzing file $id..."

#find length of alignment
length=$(bioawk -c fastx '{print length($seq)}' $faDir/$fa|head -1)

#restrict to fully aligned sites
msa_view -G ANY $faDir/$fa > $faDir/$id.fullAlign.fa
fullAlignLength=$(bioawk -c fastx '{print length($seq)}' $faDir/$id.fullAlign.fa|head -1)

#convert to nexus
python3 code/convert.py -f nexus $faDir/$id.fullAlign.fa $outNex

#find number of species
numSpecs=$(grep '>' $faDir/$id.fullAlign.fa|wc -l |awk '{print $1}')

#find proportion of missing data
propMissing=$(ruby code/get_proportion_of_missing_data.rb $outNex)
pass=$(echo $propMissing '<' .9|bc -l)


#calculate tree
#iqtree-omp -quiet -s $faDir/$id.fullAlign.fa -m GTR -pre $outDir/$id -safe -nt 4 -bb 1000 -o $outgroup
#rm -f $outDir/$id.bionj $outDir/$id.ckp.gz $outDir/$id.contree $outDir/$id.mldist $outDir/$id.splits.nex
tree=$(cat $treeDir/$id.contree)
echo $id $tree >> allTrees.tsv

# Get the proportion of missing data.
proportion_of_missing_data=`ruby code/get_proportion_of_missing_data.rb $outNex`

# Get the number of phylogenetic-informative sites.
number_of_pi_sites=`ruby code/get_number_of_pi_sites.rb $outNex`

# Get the mean bootstrap support.
bootstrap_value=`python3 code/get_mean_node_support.py $treeDir/$id.contree`
bootstrap_value_strip=`echo -n $bootstrap_value`

# Run PhiPack to test for recombination and get the resulting p-value.
p_value=`Phi -f $faDir/$fa | tail -n 2 | head -n 1 | cut -d ":" -f 2 | tr -d '[[:space:]]'`

rm -f Phi.*

# Print all output.
echo -e "$id\t$numSpecs\t$length\t$fullAlignLength\t$proportion_of_missing_data\t$number_of_pi_sites\t$bootstrap_value_strip\t$p_value" >> $statsFile

# Provide feedback on screen.
echo " done."
done
