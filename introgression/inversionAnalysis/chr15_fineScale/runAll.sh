#!/bin/bash

mkdir data
mkdir liftovers

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/windowDefinitions/HmelWindows_50KBSliding/HmelWindows.bed data

############ define variables ##############

genomes=$(halStats --genomes data/fullAlignment.hal)
startWindow=chr15_50
endWindow=chr15_84
ROIname=chr15
halFile=data/fullAlignment.hal
refGenome=HmelRef


############ run code ################
startScaf=$(awk -v win=$startWindow '$4==win' data/Hmel2Windows.bed |awk '{print $1}' )
endScaf=$(awk -v win=$endWindow '$4==win' data/Hmel2Windows.bed |awk '{print $1}' )
if [ $startScaf == $endScaf ];then
  startPos=$(awk -v win=$startWindow '$4==win' data/Hmel2Windows.bed |awk '{print $2}')
  endPos=$(awk -v win=$endWindow '$4==win' data/Hmel2Windows.bed |awk '{print $3}')
  echo $startScaf $startPos $endPos $ROIname + > $ROIname\_ROI.bed
else
  echo "Different Scaffolds"
fi

for g in $genomes
do
if [[ ! $g == *"Anc"* ]];then
code/liftoverGenome.slurm $halFile $refGenome $ROIname\_ROI.bed $g
awk '{print $4"\t"$6}' $g.$ROIname\_ROI.bed|uniq|sort|uniq|awk 'a[$1]++{ if(a[$1]==2){ print b }; print $1}; {b=$1}'|uniq > $g.candidateScafs.txt
while read line
do grep $line $g.$ROIname\_ROI.bed|
awk '{print $1"\t"$2"\t"$3"\t"$4"_"$10"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > $g.$line.bed
done < $g.candidateScafs.txt
fi
done
