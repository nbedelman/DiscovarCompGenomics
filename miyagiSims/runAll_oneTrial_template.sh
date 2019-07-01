#!/bin/bash

### take a directory with phylip files;
### generate trees and run QuIBL
### generate GENO files and get Fd
### return summary statistics that can be compiled with other replicates

####### SETUP ENVIRONMENT
mkdir GENOs
mkdir trees
mkdir logs
mkdir quibl
mkdir popStats
mkdir code

ln -s HEADDIR/code/* code/
chmod u+x code/*
export PATH=$PATH:$PWD/code

####### DEFINE VARIABLES #######
quiblDir=/n/mallet_lab/edelman/QuIBL/cython_vers/
genoDir=/n/mallet_lab/edelman/software/genomics_general

export PYTHONPATH=$PYTHONPATH:$genoDir
phymlPath=/n/sw/fasrcsw/apps/Core/phyml/2014Oct16-fasrc01/bin/phyml

phyDir=TESTDATADIR
GENODir=GENOs
combinedFile=BASENAME.geno

windSize=WINDSIZE
treeDir=trees
outgroup=4

popStatsDir=popStats
fdOut=$popStatsDir/BASENAME.popStats.csv
fdSummary=BASENAME.fdSummary.out
fdBootSummary=BASENAME.popStats.boot.summary

QUIBLDir=quibl
numLoci=5000
quiblTemplate=quiblInput.template.txt
quiblInput=$QUIBLDir/BASENAME.input.txt
treeFile=$treeDir/BASENAME.trees
quiblOut=$QUIBLDir/BASENAME.quibl.out
quiblSummary=BASENAME.quiblSummary.out
quiblBootOut=BASENAME.quiblSummary.boot.out
quiblBootSummary=BASENAME.quiblSummary.boot.summary
numBoots=100
numBootLoci=1000
baseID=BASENAME


overallOut=OVERALLOUT

###### run code ######

#convert phylips to genos and merge
convertAndMerge=`sbatch code/convertAndMerge.slurm $genoDir $phyDir $GENODir $combinedFile | cut -d " " -f 4`

getFd=`sbatch --dependency=afterok:$convertAndMerge code/getABBABABAstats.slurm $genoDir $combinedFile $windSize $fdOut $fdSummary $fdBootSummary $numBoots $numLoci|cut -d " " -f 4`

getTrees=`sbatch --dependency=afterok:$convertAndMerge code/getTrees.slurm $genoDir $phymlPath $combinedFile $treeDir $windSize|cut -d " " -f 4`
#
# #generate quibl input file
# #sed "s:TREEFILE:$treeFile:" $quiblTemplate | sed "s:OUTGROUP:$outgroup:" | sed "s:OUTPATH:$quiblOut:" > $quiblInput
#
# #runQuIBL=`sbatch --dependency=afterok:$getTrees code/runQuIBL.slurm $quiblDir $quiblInput $quiblOut $quiblSummary $numLoci|cut -d " " -f 4`
#
#
runQuIBLBoot=`sbatch --dependency=afterok:$getTrees code/bootstrapQuIBL.slurm $treeFile $quiblTemplate $outgroup $quiblDir $quiblBootOut $numBoots $numBootLoci $baseID|cut -d " " -f 4`
#
collateResults=`sbatch --dependency=afterok:$runQuIBLBoot code/collateResults.slurm $baseID $fdBootSummary $quiblBootOut $overallOut`
