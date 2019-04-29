#!/bin/bash
for filename in ./mspOut/*.txt; do
	semifile=${filename%.txt}
	stripFile=${semifile##*/}
	./seq-gen -m GTR -l 1000 -p 1000 -s 3e-9 $filename > ./sgOut/$stripFile".dat"
	back_pid=$!
	wait $back_pid
done
