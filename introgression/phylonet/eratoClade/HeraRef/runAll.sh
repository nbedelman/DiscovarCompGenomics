#!/bin/bash

##### run phyloNet_MCMC for erato clade

#going to use the already-generated fasta alignments from earlier introgression analyses.
#run 100 replicates of 100 randomly chosen alignments

##########setup environment

module purge
export MODULEPATH=$HOME/modulefiles:$MODULEPATH
module load beagle-2.1-shared
phylo=/n/mallet_lab/edelman/software/PhyloNet_3.6.1.jar



## sampling certain number of regions #########
###### define variables #######
faDir=/n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/introgression/windowTrees/eratoClade/HeraRef/10KBevery50/fastas_fullAlign/
numIntros=6
minSize=2000
maxSize=15000
numEntries=100
additionalOptions="-mc3 (2.0, 3.0) -cl 20000000"
numFiles=150
runDirectory=hotChains_150Iters_100loci_fullAlignFastas_longChain
outBase=10KBevery50.random150.fullAlign.hotChains.longChain
########## subset and format the tree file

mkdir -p $runDirectory
cp -Rf code/ $runDirectory
chmod u+x $runDirectory/code/*
cd $runDirectory
makeNexus=`sbatch code/faDirToNexus.random.new.slurm $faDir $outBase $numIntros $numFiles $numEntries $minSize $maxSize "$additionalOptions"|cut -d " " -f 4`
let fileNameNum="$numFiles - 1"
for i in $(seq 0 $fileNameNum)
do outDir=$outBase.s$i\_outDir
mkdir $outDir
sbatch --dependency=afterok:$makeNexus code/runPhyloNet.slurm $phylo $outBase.s$i.nexus $outBase.s$i.out
done
