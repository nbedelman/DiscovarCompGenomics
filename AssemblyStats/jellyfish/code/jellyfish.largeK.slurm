#!/bin/bash

#SBATCH -J jellyfish
#SBATCH -n 32
#SBATCH -t 0-05:00
#SBATCH --mem=100000
#SBATCH -p serial_requeue
#SBATCH -e jellyfish_%A.err
#SBATCH -o jellyfish_%A.out

reads1=$1
reads2=$2
outBase=$3

jellyfish count -m 31 -o $outBase.out -c 3 -s 100M -t 32 -F 2 <(zcat $reads1) <(zcat $reads2)
jellyfish histo -t 32 $outBase.out > $outBase.hist
