#!/bin/bash

###################### Setup Environment #################

mkdir -p data
mkdir -p code
mkdir -p mafs
mkdir -p fastas
mkdir -p VCFs

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal
ln -s /n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Hmel2.fa_1kb.fa data/refGenome.fa
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/Hmel2pt5/Hmel2.transitions.tsv data/mapFile.tsv

chmod u+x code/*


module load gcc/6.3.0-fasrc01 AdmixTools/20171018-fasrc01

################ Define Variables ###############
refGenome=HmelRef
chromMap=data/mapFile.tsv

################# Run Code ######################
#first, extract MAF for each chromosome
chroms=$(grep ">" data/refGenome.fa |cut -d ">" -f 2|cut -d " " -f 1|tail -573)
for chrom in $chroms
do extract=`sbatch code/extractChromMaf.slurm data/fullAlignment.hal $refGenome $chrom | cut -d " " -f 4`
echo $extract: >> extractJobs.txt
mv $chrom\_to$refGenome.maf mafs
done

#########re-ran, restricting analysis to only single copy. did:
for maf in mafs
do code/ $maf singleCopyMafs/$(basename $maf .maf)
done
#and then ran the rest of the code on

#then, convert maf to fasta
extractJobs=$(tr -d "\n" < extractJobs.txt)
extractJobs=$(basename $extractJobs :)
maf2Fasta=`sbatch --dependency=afterok:$extractJobs code/mafToFasta.slurm mafs| cut -d " " -f 4`

#then, convert FASTA to VCF and Eigenstrat
fasta2Eigenstrat=`sbatch --dependency=afterok:$maf2Fasta code/fastaToEigenstrat.slurm fastas| cut -d " " -f 4`

#In this case, since we're mapping to Hmel and have a reference chromosomal order, convert the scaffold-level eigenstrat files to chromosomes.
# code/scafNamesToNums.py allGenomes.snp allGenomes.geno allHeliconiini_toHmel_placed $chromMap
#
# qpDstat -p parFiles/parDstat_allHeliconius > allHeliconius.D.out
