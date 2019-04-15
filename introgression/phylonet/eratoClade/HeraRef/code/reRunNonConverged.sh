#!/bin/bash

#after a run of PhyloNet is complete, re-start any iterations that do not have a single network with >= 75% PP.

module purge
#export MODULEPATH=$HOME/modulefiles:$MODULEPATH
module load centos6/0.0.1-fasrc01 intel/13.0.079-fasrc01 beagle/2.1.trunk-fasrc01 jdk
phylo=/n/mallet_lab/edelman/software/PhyloNet_3.6.7.jar

for i in *.s*.out
do if grep -q "Rank = 0" $i;then
prob=$(grep "Rank = 0" $i|cut -d ";" -f 3|cut -d "=" -f 2|cut -d . -f 1)
if [ $prob -le 74 ];then
code/makeResumingNexus.py $(basename $i .out).nexus $(basename $i .out)_outDir $(basename $i .out).resume.nexus
sbatch code/runPhyloNet.slurm $phylo $(basename $i .out).resume.nexus $(basename $i .out).resume.out
mv $i $i.old
fi
else
  code/makeResumingNexus.py $(basename $i .out).nexus $(basename $i .out)_outDir $(basename $i .out).resume.nexus
  sbatch code/runPhyloNet.slurm $phylo $(basename $i .out).resume.nexus $(basename $i .out).resume.out
  mv $i $i.old
fi
done
