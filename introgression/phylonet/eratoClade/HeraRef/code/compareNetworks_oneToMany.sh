#!/bin/sh

#run the phylonet cmpnets command to get a pairwise metric of similarity between one tree and each of a list of trees.

module purge
module load java/1.8.0_45-fasrc01
module load beagle/2.1.trunk-fasrc03
phylo=/n/mallet_lab/edelman/software/PhyloNet_3.6.1.jar

qpNetwork=firstPhyloTree.nwk
phyloNetworks=bestTrees.tree
nexBase=firstPhyloToPhyloComp
calc=tri

code/makeCompNexus.py $qpNetwork $phyloNetworks $nexBase $calc
for i in "$nexBase"*.nexus
do java -jar $phylo $i >> $nexBase.allComparisons.$calc.out
done
