#!/bin/bash

#SBATCH -J filterAlignments
#SBATCH -p serial_requeue,shared,general
#SBATCH -t 1-00:00
#SBATCH --mem=5000
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -e filterAlignments.err
#SBATCH -o filterAlignments.out

faDir=$1
outDir=$2
statsFile=$3
outgroup=$4


mkdir -p $outDir

#Calculate trees and output relevant stats
for fa in $(ls $faDir)
do

id=$(basename $fa .fa)
outNex=$outDir/$id.nex

# Provide feedback on screen.
#echo -n "Analyzing file $id..."

#convert to nexus
python3 code/convert.py -f nexus $faDir/$fa $outNex

#find length of alignment
length=$(bioawk -c fastx '{print length($seq)}' $faDir/$fa|head -1)

#find number of species
numSpecs=$(grep '>' $faDir/$fa|wc -l |awk '{print $1}')

#find proportion of missing data
propMissing=$(ruby code/get_proportion_of_missing_data.rb $outNex)
pass=$(echo $propMissing '<' .9|bc -l)
#calculate tree
iqtree-omp -quiet -s $faDir/$fa -m GTR -pre $outDir/$id -safe -nt 4 -bb 1000 -o $outgroup
rm -f $outDir/$id.bionj $outDir/$id.ckp.gz $outDir/$id.contree $outDir/$id.mldist $outDir/$id.splits.nex
tree=$(cat $outDir/$id.treefile)
echo $id $tree >> allTrees.tsv

# Get the proportion of missing data.
proportion_of_missing_data=`ruby code/get_proportion_of_missing_data.rb $outNex`

# Get the number of phylogenetic-informative sites.
number_of_pi_sites=`ruby code/get_number_of_pi_sites.rb $outNex`

# Get the mean bootstrap support.
bootstrap_value=`python3 code/get_mean_node_support.py $outDir/$id.treefile`
bootstrap_value_strip=`echo -n $bootstrap_value`

# Run PhiPack to test for recombination and get the resulting p-value.
p_value=`Phi -f $faDir/$fa | tail -n 2 | head -n 1 | cut -d ":" -f 2 | tr -d '[[:space:]]'`

rm -f Phi.*

# Print all output.
echo -e "$id\t$numSpecs\t$length\t$proportion_of_missing_data\t$number_of_pi_sites\t$bootstrap_value_strip\t$p_value" >> $statsFile

# Provide feedback on screen.
echo " done."
done
