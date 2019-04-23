#!/bin/bash

fastqDir=$1

module load jellyfish

for f in $(ls $fastqDir/*.fastq.gz|cut -d "." -f 1|uniq)
do
sbatch code/jellyfish.slurm $f.1.fastq.gz $f.2.fastq.gz $(basename $f)
#code/jellyfish.slurm $f.1.fastq.gz $f.2.fastq.gz $(basename $f)
done &
