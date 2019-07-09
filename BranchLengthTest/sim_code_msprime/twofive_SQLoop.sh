#!/bin/bash
for folder in ./*/250000; do
	upper=$(basename $folder)
	for filename in $folder/*; do
		#echo ./sgRepeat/${folder#./}/$(basename $filename)
		mkdir -p ./sgRepeat/${folder#./}/$(basename $filename)
		for file in $filename/*; do
			#echo ./sgRepeat/${folder#./}/$(basename $filename)/$(basename $file .txt).dat
			./seq-gen -m GTR -l 1000 -p 1 -s 3e-9 $file > ./sgRepeat/${folder#./}/$(basename $filename)/$(basename $file .txt).dat
		done

		#mkdir -p ./sgOut/$upper/$(basename $subfold)/
		#../seq-gen -m GTR -l 1000 -p 1 -s 3e-9 ./mspOut/$upper/$(basename $subfold)/$(basename $filename) > ./sgOut/$upper/$(basename $subfold)/$(basename $filename .txt).dat
			
	done
		#echo $upper/$(basename $subfold)/
		#../seq-gen -m GTR -l 1000 -p 1 -s 3e-9 $
done
	#mkdir ./sgOut/$folder
	

#for filename in ./mspOut/*.txt; do
#	semifile=${filename%.txt}
#	stripFile=${semifile##*/}
#	./seq-gen -m GTR -l 1000 -p 1 -s 3e-9 $filename > ./sgOut/$stripFile".dat"
#	back_pid=$!
#	wait $back_pid
#done
