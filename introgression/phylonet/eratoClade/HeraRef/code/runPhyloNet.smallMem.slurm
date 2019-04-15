#!/bin/bash

#SBATCH -J PhyloNet
#SBATCH -p general,shared
#SBATCH -n 32
#SBATCH --mem=50000
#SBATCH -t 0-01:00
#SBATCH -e runPhyloNet.err
#SBATCH -o runPhyloNet.out


ulimit -c unlimited

phyloFile=$1
nexusFile=$2
outFile=$3

java -jar $phyloFile $nexusFile > $outFile
