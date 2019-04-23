#!/bin/bash

####Get final stats on all alignments/trees actually used to generate species trees

############################################
################ Setup Environment #########
############################################

chmod u+x code/*
export PATH=code/:$PATH

module load bedtools2
module load phast
module load Anaconda
module load python/3.6.3-fasrc02
module load ucsc
module load centos6/0.0.1-fasrc01
module load R/3.4.1-fasrc01

mkdir -pv ~/apps/R
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER

############################################
################ Define Variables ##########
############################################

allTreesDir=/n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/phylogenies/phylogeny_trees_fastas_FINAL/single_locus_trees
allFastasDir=/n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/phylogenies/phylogeny_trees_fastas_FINAL/fasta_alignments
statsBase=singleLocusStats

############################################
################ Run Code ##################
############################################


#### make a directory for each subset to compute basic statistics
for treeDir in $allTreesDir/*
do
  blockSet=$(basename $treeDir)
  mkdir -p $blockSet
  cp -R code $blockSet
  cd $blockSet
  statsFile=$statsBase.$blockSet.txt
  echo -e "block_id\tnumSpecies\tfullLength\tfullyAlignedSites\tp_missing\tn_pi_sites\tmean_bootstrap_support\tphi_p" > $statsFile
  total=$(ls $treeDir|wc -l)
  start=0
  while [ $start -le $total ];
    do
    faDir=fastas_$start
    trees=trees_$start
    outDir=output_$start
    mkdir $faDir
    mkdir $trees
    mkdir $outDir
    theseFastas=$(ls $allFastasDir/$blockSet|tail -n +$start|head -499)
    for file in $theseFastas;
      do
      fileBase=$(basename $file .fas)
        ln -s $allFastasDir/$blockSet/$fileBase.fas $faDir
        ln -s $allTreesDir/$blockSet/$fileBase.contree $trees
      done
    getStats=`sbatch code/filterAlignments.slurm $faDir $statsFile $trees $outDir|cut -d " " -f 4`
    let start="$start + 500"
  done
  R CMD BATCH --no-save --no-restore "--args $statsFile $blockSet" code/blockStatisticVisualization.R
  cd ..
done
