#!/bin/bash

#SBATCH -J MS_martinSims
#SBATCH -p serial_requeue,shared,general
#SBATCH -t 0-03:00
#SBATCH --mem=1000
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -e MS_martinSims.err
#SBATCH -o MS_martinSims.out

#Just going to run Simon's simulations to make sure I'm not doing anything stupid with fd and D calcs


######setup directory ######
mkdir msTrees
mkdir seqs
mkdir genos
mkdir infTreeDir
mkdir quibl
mkdir logs

module load  centos6/0.0.1-fasrc01 ms/2014-fasrc01
module load phyml/2014Oct16-fasrc01
module load Anaconda/5.0.1-fasrc01
source activate quibl

###### define variables #######
genoDir=/n/mallet_lab/edelman/software/genomics_general
phymlPath=/n/sw/fasrcsw/apps/Core/phyml/2014Oct16-fasrc01/bin/phyml
quiblPath=/n/mallet_lab/edelman/QuIBL/cython_vers/
export PYTHONPATH=$PYTHONPATH:$genoDir

numLoci=1000
msTreeDir=msTrees
seqDir=seqs
nonIntroFrac=NONINTROFRAC
introTime=INTROTIME
base=BASENAME

GENODir=genos
combinedFile=$base.combined.geno
inferredTreeDir=infTreeDir
windSize=5000
quiblTemplate=/n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/miyagiSims/msSims/quiblInput.template.txt

header=$(echo \#CHROM POS $(seq 1 4 | tr "\n" " "))
echo $header > $combinedFile

#######

#code/simulateSeqs.slurm $numLoci $msTreeDir $seqDir $base $nonIntroFrac $introTime

#code/runPopStats.slurm $seqDir $genoDir $GENODir $windSize $base $combinedFile

code/makeTreesAndQuibl.sh $genoDir $combinedFile $inferredTreeDir $phymlPath $quiblTemplate $quiblPath
