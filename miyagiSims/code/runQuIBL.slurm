#!/bin/bash

#SBATCH -J runQuIBL
#SBATCH -n 1
#SBATCH -p serial_requeue,shared,general
#SBATCH -t 0-10:00
#SBATCH --mem=5000
#SBATCH -e logs/runQuIBL.err
#SBATCH -o logs/runQuIBL.out

quiblLoc=$1
inputFile=$2
quiblOut=$3
quiblSummary=$4
numLoci=$5

python $quiblLoc/QuIBL_cyth.py $inputFile
sed 's/,/\t/g' $quiblOut |awk '$1=="1_2_3"||$1=="1_3_2"|| $1=="2_3_1"||$1=="2_1_3"||$1=="3_2_1"||$1=="3_1_2"'|awk -v nL=$numLoci '$2=="1" {print $4","$7","$6*$11/nL","$10-$9}' >> $quiblSummary
