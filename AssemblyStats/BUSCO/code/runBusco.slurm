#!/bin/bash

#SBATCH -J runBusco
#SBATCH -p serial_requeue
#SBATCH -n 2
#SBATCH -c 16
#SBATCH -t 0-02:00
#SBATCH --mem=50000
#SBATCH -e runBusco%A.err
#SBATCH -o runBusco%A.out

mkdir -p allSummaries

genome=$1
lineage=$2

BUSCO.py -i $genome -c 32 -m geno -l $lineage -sp heliconius_melpomene1  -o $(basename $genome .fa) -r -f
cp run_$(basename $genome .fa)/short_summary* allSummaries
summary=$(grep "C:" run_$(basename $genome .fa)/short_summary*)
echo $genome $summary >> allSummaries.txt
