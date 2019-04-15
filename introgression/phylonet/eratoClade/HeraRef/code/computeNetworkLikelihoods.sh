#!/bin/bash
#takes a list of networks (i.e. best network from MCMC output)
#and a list of gene trees
#creates a newick file to run PhyloNet's "CalGTProb", and runs it.
#returns a file with the networks along with their likelihood values.

######### SETUP ENVIRONMENT ######
#don't need to do this if it's running as part of a pipeline that already did it
module purge
module load java/1.8.0_45-fasrc01
module load beagle/2.1.trunk-fasrc03
phylo=/n/mallet_lab/edelman/software/PhyloNet_3.6.1.jar

########## DEFINE VARIABLES ########
geneTrees=/n/regal/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/introgression/eratoClade_5KBEvery50/phyMLTrees.ge1KB.rooted.newickOnly.txt
networks=hotChainTrial.tree
summaryFile=hotChainTrial.networkLikelihoods.out
########## RUN SCRIPTS #######
num=1
while read line
do nexFile=hotChainTrial.network$num.nexus
outFile=hotChainTrial.network$num.out
code/makeCalGTProbNexus.py $nexFile $line $geneTrees $additionalOptions
sbatch code/runPhyloNet.smallMem.slurm $phylo $nexFile $outFile
let num="$num+1"
done < $networks
