#!/bin/bash

#SBATCH -J fastaToEigenstrat
#SBATCH -p general,shared,serial_requeue
#SBATCH -t 0-01:00
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -e fastaToEigenstrat.err
#SBATCH -o fastaToEigenstrat.out

module load Anaconda

source activate my_root

export PYTHONPATH=$PYTHONPATH:/n/mallet_lab/edelman/software/gdc
#
module load phast/1.4-fasrc01
module load vcftools

fastaDir=$1
base=$2

mkdir -p $base\_VCFs

for fastaFile in $fastaDir/*.fa
do numSeqs=$(grep ">" $fastaFile|wc -l)
if [ $numSeqs -ge 17 ]; then
VCFfile=$base\_VCFs/$(basename $fastaFile .fa).vcf
scaffold=$(basename $fastaFile|cut -d "_" -f 1|cut -d "." -f 1)
msa_view --seqs Etal,Ldor,Htel,Hhsa,HeraRef,HeraDisco,Hhim,HeraHhimHyb,Hsar,Hdem,Hbes,Hpar,Hnum,Htim,Hcyd,HmelRef,HmelDisco $fastaFile|
sed 's/ //g' > $fastaFile\_tmp
numSeqs=$(grep ">" $fastaFile\_tmp|wc -l)
if [ $numSeqs -eq 17 ]; then
snp-sites -cv -o $VCFfile $fastaFile\_tmp
sed -i "s/^1/$scaffold/g" $VCFfile

###This part is specific to Hmel2###
chrom=${scaffold:0:7}
awk -v chr=$chrom '$1=chr' $VCFfile|awk $3=$1:$2
####################################
fi
fi
rm -f $fastaFile\_tmp
done

vcf-concat $base\_VCFs/* > allGenomes.$base.vcf

python code/vcfToEigenstrat.py -v allGenomes.$base.vcf -o allGenomes.$base
