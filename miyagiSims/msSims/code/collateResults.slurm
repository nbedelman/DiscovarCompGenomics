#!/bin/bash

#SBATCH -J collateResults
#SBATCH -p serial_requeue,shared,general
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 0-00:10
#SBATCH --mem=100
#SBATCH -e logs/collateResults.err
#SBATCH -o logs/collateResults.out

id=$1
fdResults=$2
quiblResults=$3
totalOut=$4

fd=$(head -1 $fdResults)

python code/getBootSummaryStats.py $quiblResults $(basename $quiblResults .out).summary

quibl=$(head -1 $(basename $quiblResults .out).summary )

IDName=$(echo $id|rev|cut --complement -d _ -f-2|rev)
IDInfo=$(echo $id|rev|cut -d _ -f-2|rev|sed 's/_/,/g')

echo $IDName,$IDInfo,$fd,$quibl >> $totalOut
