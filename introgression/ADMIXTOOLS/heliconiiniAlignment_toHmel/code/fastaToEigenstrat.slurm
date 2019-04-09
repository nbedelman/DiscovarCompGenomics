#!/bin/bash

#SBATCH -J fastaToEigenstrat
#SBATCH -p general
#SBATCH -t 0-01:00
#SBATCH -n 1
#SBATCH --mem=5000
#SBATCH -e fastaToEigenstrat.err
#SBATCH -o fastaToEigenstrat.out

module load Anaconda

source activate my_root

export PYTHONPATH=$PYTHONPATH:/n/mallet_lab/edelman/software/gdc
#
module load gcc/4.8.2-fasrc01 phast
module load vcftools

fastaDir=$1

for fastaFile in $fastaDir/*.fa
do numSeqs=$(grep ">" $fastaFile|wc -l)
if [ $numSeqs -ge 22 ]; then
VCFfile=VCFs/$(basename $fastaFile .fa).vcf
scaffold=$(basename $fastaFile|cut -d "_" -f 1)
msa_view --seqs  Eueides_tales_helico3,Heliconius_telesiphe_helico3,Heliconius_hecalesia_helico3,Heliconius_himera_helico3,Heliconius_himera_father_helico3,Heliconius_erato_x_himera_F1_helico3,Heliconius_erato_demophoon_v1,Heliconius_erato_lativitta_v3,Heliconius_erato_mother_helico3,Heliconius_demeter_helico3,Heliconius_sara_helico3,Laparus_doris_helico3,Heliconius_burneyi_helico3,Heliconius_melpomene_helico3,Heliconius_melpomene_melpomene_hmel2,Heliconius_cydno_helico3,Heliconius_timareta_helico3,Heliconius_numata_helico3,Heliconius_besckei_helico3,Heliconius_hecale_helico3,Heliconius_elevatus_helico3,Heliconius_pardalinus_helico3 $fastaFile|
sed 's/ //g' > $fastaFile\_tmp
numSeqs=$(grep ">" $fastaFile\_tmp|wc -l)
if [ $numSeqs -eq 22 ]; then
snp-sites -cv -o $VCFfile $fastaFile\_tmp
sed -i "s/^1/$scaffold/g" $VCFfile
fi
fi
rm -f $fastaFile\_tmp
done

vcf-concat VCFs/* > allGenomes.vcf

python code/vcfToEigenstrat.py -v allGenomes.vcf -o allGenomes
