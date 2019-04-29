#!/bin/bash
COUNTER=$1
for filename in ./mspOut/*.txt; do
	semifile=${filename%.txt}
	stripFile=${semifile##*/}
	./seq-gen -m GTR -l 1000 -p 1000 -s 3e-9 $filename > ./sgOut/$stripFile".dat"
	back_pid=$!
	wait $back_pid
done
for filename in ./sgOut/*.dat; do
	semifile=${filename%.dat}
	stripFile=${semifile##*/}
	./PhyML3164 -q -i $filename -s BEST -o 'tlr' > ./pmlOut/$stripFile
	#./PhyML3164 -q -i $filename > ./pmlOut/$stripFile
	back_pid=$!
	wait $back_pid
done
cat ./sgOut/*tree* > ./fullTrees/treeCat.$COUNTER.txt
python ./rootingGTtest.py ./fullTrees/treeCat.$COUNTER.txt ./tripleSets/triple$COUNTER
python ~/Documents/GT_Wakeley/FIXED_LAMBDA_PLexMax_Cythonize/cyth_FL_exMax.py ./tripleSets/triple$COUNTER.p ./outFiles/out$COUNTER 2 20 0.006
rm ./sgOut/*
rm ./pmlOut/*
