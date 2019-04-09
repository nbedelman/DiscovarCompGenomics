#!/bin/bash

#SBATCH -J qpGraph
#SBATCH -p general
#SBATCH -t 0-01:00
#SBATCH -n 1
#SBATCH --mem=50000
#SBATCH -e qpGraph.err
#SBATCH -o qpGraph.out

#module load gcc/6.3.0-fasrc01 AdmixTools/20171018-fasrc01

parFile=$1
graph=$2
outBase=$3

qpGraph -p $parFile -g $graph -o $outBase.outGraph -d $outBase.dot > $outBase.out
dot -Tps $outBase.dot > $outBase.dot.ps
