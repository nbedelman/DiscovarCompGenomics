#!/bin/bash

#SBATCH -J runqpDstat
#SBATCH -p serial_requeue,general,shared
#SBATCH -n 1
#SBATCH -t 0-01:00
#SBATCH --mem=50000
#SBATCH -e Dstat.err
#SBATCH -o Dstat.out

module load gsl
module load CLAPACK/3.2.1-fasrc01
module load intel/17.0.4-fasrc01 OpenBLAS/0.2.20-fasrc01
export PATH=$PATH:/n/mallet_lab/edelman/software/AdmixTools/bin

parFile=$1
outFile=$2

qpDstat -p $parFile > $outFile
