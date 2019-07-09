#!/bin/bash
for folder in ./mspOut/*/; do
		upper=$(basename $folder)
	for subfold in $folder*/; do
		for filename in $subfold/*; do
			mkdir -p ./sgOut/$upper/$(basename $subfold)/
			../seq-gen -m GTR -l 5000 -p 1 -s 3e-9 ./mspOut/$upper/$(basename $subfold)/$(basename $filename) > ./sgOut/$upper/$(basename $subfold)/$(basename $filename .txt).dat
			
		done
		#echo $upper/$(basename $subfold)/
		#../seq-gen -m GTR -l 1000 -p 1 -s 3e-9 $
	done
done
	#mkdir ./sgOut/$folder
	

#for filename in ./mspOut/*.txt; do
#	semifile=${filename%.txt}
#	stripFile=${semifile##*/}
#	./seq-gen -m GTR -l 1000 -p 1 -s 3e-9 $filename > ./sgOut/$stripFile".dat"
#	back_pid=$!
#	wait $back_pid
#done
