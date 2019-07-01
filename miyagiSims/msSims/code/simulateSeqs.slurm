#!/bin/bash

numLoci=$1
treeDir=$2
seqDir=$3
base=$4
nonIntro=$5
introTime=$6

for rep in $(seq 1 $numLoci)
do
  treeFile=$treeDir/$base.$rep.tree
  seqFile=$seqDir/$base.$rep.phy

  ms 4 1 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -es $introTime 2 $nonIntro -ej $introTime 5 3 -r 50 5000 -T | tail -n +4 | grep -v // > $treeFile
  partitions=($(wc -l $treeFile))

  seq-gen -mHKY -l 5000 -s 0.01 -p $partitions < $treeFile > $seqFile

done


# ms 32 1 \ #32 samples, 1 chrom each
# -I 4 8 8 8 8 \ #4 populations of 8 individuals each
# -ej 1 2 1 \ #at time 1, combine pop2 into pop1
# -ej 2 3 1 \ #at time 2, combine pop3 into pop1
# -ej 3 4 1 \ #at time 3, combine pop4 into pop1
# -es 0.1 2 0.8 \ #at time 0.1, split pop 2: 80% stay in pop2, 20% move into pop5
# -ej 0.1 5 3 \ #also at time 0.1,move all those pop 5 individuals into pop3. Forward in time, these moves correspond to 20% migration of 3 into 2
# -r 50 5000 \ #rho=50. for Ne=1e6 and seq length 5KB, this means per-bp r=1x10^-8
# -T | \ #return trees
