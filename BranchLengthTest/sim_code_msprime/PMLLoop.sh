#!/bin/bash
for filename in ./sgOut/*.dat; do
	semifile=${filename%.dat}
	stripFile=${semifile##*/}
	./PhyML3164 -q -i $filename -s BEST -o 'tlr' > ./pmlOut/$stripFile
	#./PhyML3164 -q -i $filename > ./pmlOut/$stripFile
	back_pid=$!
	wait $back_pid
done
