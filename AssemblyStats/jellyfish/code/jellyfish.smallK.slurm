#!/bin/bash

#SBATCH -J jellyfish
#SBATCH -n 2
#SBATCH -c 16
#SBATCH -t 0-01:00
#SBATCH --mem=500000
#SBATCH -p serial_requeue
#SBATCH -e jellyfish_%A.err
#SBATCH -o jellyfish_%A.out

reads1=$1
#reads2=$2
outBase=$2

jellyfish count -m 15 -o $outBase.out -c 3 -s 100M -t 32 -F 2 <(zcat $reads1) #<(zcat $reads2)
cat $outBase.out* > $outBase.out
jellyfish histo -t 32 $outBase.out > $outBase.hist
