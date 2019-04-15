#!/bin/bash

#template to generate trees from windows

############################################
################ Setup Environment #########
############################################

mkdir -p data
mkdir -p code
mkdir -p mafs
mkdir -p fastas
mkdir -p beds
mkdir -p trees

chmod u+x code/*
export PATH=code/:$PATH

module load bedtools2
module load phast
module load Anaconda
module load python/3.6.3-fasrc02
module load centos6/0.0.1-fasrc01  gcc/7.1.0-fasrc01 iqtree/1.5.5-fasrc01
module load ucsc

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal

##Change this depending on outgroup and window size
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/windowDefinitions/HeraWindows_10KBAbutting/HeratoWindows.bed data/windows.bed


############################################
################ Define Variables ##########
############################################

#extract maf and fasta alignments for each window
refGenome=HeraRef
halPath=data/fullAlignment.hal
targetGenomes=HmelRef,HeraRef,Hhim,Hhsa,Htel,Hdem,Hsar
windowBed=data/windows.bed

#make and summarize trees
faDir=fastas
treeDir=$(pwd)/trees
statsFile=heliconiini_HeraRef_windowStats.txt
alignDir=output
outgroup=HmelRef

############################################
################ Run Code ##################
############################################

chroms=$(awk '{print $4}' $windowBed|cut -d "_" -f 1 |sort|uniq)
for chrom in $chroms
do grep $chrom\_ $windowBed > $refGenome.$chrom.windows.bed
thisJob=`sbatch code/halToMaf_regions.slurm $refGenome $halPath $targetGenomes $refGenome.$chrom.windows.bed | cut -d " " -f 4`
echo $thisJob: >> halToMaf.jobs.txt
done

halToMafJobs=$(tr -d "\n" < halToMaf.jobs.txt)
halToMafJobs=$(basename $halToMafJobs :)

#setup stats file
echo -e "block_id\tnumSpecies\tfullLength\tp_missing\tn_pi_sites\tmean_bootstrap_support\tphi_p" > $statsFile

for chrom in {1..21}
do getTrees=`sbatch --dependency=afterok:$halToMafJobs code/filterAlignments.slurm $chrom $faDir chr$chrom\_$alignDir $statsFile $outgroup|cut -d " " -f 4`
echo $getTrees:>>treeJobs.txt
done
treeJobs=$(tr -d "\n" < treeJobs.txt)
treeJobs=$(basename $treeJobs :)
