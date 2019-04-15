#!/bin/bash

#makeReferenceWindows.slurm


#SBATCH -n 1
#SBATCH -t 0-10:00
#SBATCH --mem=1000
#SBATCH -p general,serial_requeue,shared
#SBATCH -o makeReferenceWindows.out
#SBATCH -e makeReferenceWindows.err


#This really only works for Hmel2.5-Hmel2 because the TSV is not a standard format.
transitionTSV=$1
transitionBase=$(basename $transitionTSV .tsv)
newRefBed=$2
newRefBase=$(basename $newRefBed .bed)
windowSize=$3
slide=$4

#create windows for the downstream analysis
chroms=$(awk '{print $1}' $transitionTSV |sort|uniq)
for chrom in $chroms
do
largest=$(awk -v chrom=$chrom '$1==chrom {print}' $transitionTSV|sort -nrk3,3 |head -1 |awk '{print $3}' )
echo chr$chrom 1 $largest chr$chrom >> $transitionBase\_fullChroms.bed
done
sed -e 's/ [ ]*/\t/g' $transitionBase\_fullChroms.bed > $transitionBase\_fullChroms_tab.bed
sort -k1 $transitionBase\_fullChroms_tab.bed > $transitionBase\_fullChroms.bed && rm -f $transitionBase\_fullChroms_tab.bed
bedtools makewindows -b $transitionBase\_fullChroms.bed -w $windowSize -s $slide -i srcwinnum> $transitionBase\_Windows.bed

awk '{print "chr"$1"\t"$2"\t"$3"\t"$7"\t.\t"$10}' $transitionTSV > $transitionBase.bed
bedtools intersect -a $transitionBase\_Windows.bed -b $transitionBase.bed -wao > $newRefBase.tmp

windows=$(awk '{print $4}' $transitionBase\_Windows.bed)
for window in $windows
do
numSeqs=$(awk -v window=$window '$4==window' $newRefBase.tmp|wc -l)
if [ $numSeqs == 1 ];then
  #in this case, the full 50KB window is within a single scaffold. just have to figure out which 50KB section it is.
  awk -v window=$window -v winSize=$windowSize '$4==window {print $8"\t"$2-$6+1"\t"$2-$6+winSize+1"\t"$4"\t.\t"$10}' $newRefBase.tmp >> $newRefBed
else
first=$(awk -v window=$window '$4==window' $newRefBase.tmp |sort -nk 6,6|head -1)
firstLength=$(echo $first|awk '{print $7-$6}')
firstOverlap=$(echo $first|awk '{print $11}')
if [ $firstLength == $firstOverlap ];then
  echo $first | awk '{print $8"\t"1"\t"$7-$6+1"\t"$4"\t.\t"$10}' >> $newRefBed
else
  echo $first| awk -v overlap=$firstOverlap '{print $8"\t"$7-$6+1-overlap"\t"$7-$6+1"\t"$4"\t.\t"$10}' >> $newRefBed
fi

if [ $numSeqs -ge 2 ];then
  awk -v window=$window '$4==window' $newRefBase.tmp|sort -nk 6,6|head -n -1|tail -n +2 |awk '{print $8"\t"1"\t"$7-$6+1"\t"$4"\t.\t"$10}'>> $newRefBed
fi

last=$(awk -v window=$window '$4==window' $newRefBase.tmp |sort -nk 6,6|tail -1)
lastLength=$(echo $last|awk '{print $7-$6}')
lastOverlap=$(echo $last|awk '{print $11}')
if [ $lastLength == $lastOverlap ];then
  echo $last | awk '{print $8"\t"1"\t"$7-$6+1"\t"$4"\t.\t"$10}' >> $newRefBed
else
  echo $last | awk -v overlap=$lastOverlap '{print $8"\t"1"\t"1+overlap"\t"$4"\t.\t"$10}' >> $newRefBed
fi
fi
done
