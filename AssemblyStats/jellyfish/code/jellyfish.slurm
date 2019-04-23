#!/bin/bash

#SBATCH -J jellyfish
#SBATCH -n 4
#SBATCH -c 4
#SBATCH -t 0-05:00
#SBATCH --mem=100000
#SBATCH -p general,shared,serial_requeue
#SBATCH -e jellyfish_%A.err
#SBATCH -o jellyfish_%A.out

reads1=$1
reads2=$2
outBase=$3

jellyfish count -m 15 -o $outBase.out -s 2000M -t 16 -F 2 -C <(zcat $reads1) <(zcat $reads2)
jellyfish histo -t 16 $outBase.out > $outBase.hist
