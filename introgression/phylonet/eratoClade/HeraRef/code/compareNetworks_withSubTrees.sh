#!/bin/sh

#run the phylonet charnet command to get an output of the subtrees in each network.
#Then, compare those subtrees to get an idea about which evolutionary paths are most well supported.

module purge
module load Anaconda
module load java/1.8.0_45-fasrc01
module load beagle/2.1.trunk-fasrc03
phylo=/n/mallet_lab/edelman/software/PhyloNet_3.6.1.jar
source activate my_root

qpNetwork=qpGraphTree.nwk
phyloNetworks=fastas.10KBevery50.random1000.credibleTrees.noWeights.out
nexBase=subTreeComp
calc=tree

code/makeCharNexus.py $qpNetwork $nexBase $calc
code/makeCharNexus.py $phyloNetworks $nexBase $calc
for i in "$nexBase"*.nexus
do echo $i >> $nexBase.allNets.$calc.out
java -jar $phylo $i |tail -n +3>> $nexBase.allNets.$calc.out
done

code/compareSubTrees.py $nexBase.allNets.$calc.out $nexBase.allNets.$calc
