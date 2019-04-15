#!/bin/bash

#SBATCH -J PhyloNet
#SBATCH -p serial_requeue
#SBATCH -n 2
#SBATCH -c 16
#SBATCH --mem=50000
#SBATCH -t 1-00:00
#SBATCH -e runPhyloNet.err
#SBATCH -o runPhyloNet.out


# ulimit -c unlimited

phyloFile=$1
nexusFile=$2
outFile=$3

java -jar -Xmx10G $phyloFile $nexusFile > $outFile
