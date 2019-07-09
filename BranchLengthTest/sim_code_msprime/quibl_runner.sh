#!/bin/bash
for folder in ./mspOut/*/; do
		upper=$(basename $folder)
	for subfold in $folder*/; do
		for filename in $subfold/*; do
			mkdir -p ./quiblOut/$upper/$(basename $subfold)/
#			../seq-gen -m GTR -l 1000 -p 1 -s 3e-9 ./mspOut/$upper/$(basename $subfold)/$(basename $filename) > ./sgOut/$upper/$(basename $subfold)/$(basename $filename .txt).dat
			
		done
		#echo $upper/$(basename $subfold)/
		#../seq-gen -m GTR -l 1000 -p 1 -s 3e-9 $
	done
done
