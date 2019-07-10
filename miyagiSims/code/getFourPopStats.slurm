#!/bin/bash

#SBATCH -J getFourPopStats
#SBATCH -p serial_requeue,shared,general
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 0-01:00
#SBATCH --mem=1000
#SBATCH -e logs/getFourPopStats.err
#SBATCH -o logs/getFourPopStats.out

genoDir=$1
genoFile=$2
windSize=$3
outFile=$4

python $genoDir/fourPopWindows.py -w $windSize -g $genoFile -f phased -P1 1 1 -P2 2 2 -P3 3 3 -O 4 4 -T 16 -o $outFile
