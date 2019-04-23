#!/bin/bash

mkdir data
mkdir liftovers

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/windowDefinitions/HmelWindows_50KBSliding/HmelWindows.bed  data

############ define variables ##############

genomes=$(halStats --genomes data/fullAlignment.hal)
startWindow=chr2_170
endWindow=chr2_260
ROIname=chr2
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

for g in HmelRef HmelDisco HeraRef HeraDisco Hdem Hsar Hhsa Htel Hhim Etal #$genomes
do
if [[ ! $g == *"Anc"* ]];then
code/liftoverGenome.slurm $halFile $refGenome $ROIname\_ROI.bed $g
awk '{print $4"\t"$6}' $g.$ROIname\_ROI.bed|uniq|sort|uniq|awk 'a[$1]++{ if(a[$1]==2){ print b }; print $1}; {b=$1}'|uniq > $g.candidateScafs.txt

while read line
# do grep $line $g.$ROIname\_ROI.bed > $g.$line.bed
do hal2maf $halFile $g.$line.maf --refGenome $g --refSequence $line --targetGenomes $g,$refGenome
code/mafToBed.py $g.$line.maf $g.$line.bed
grep Hmel $g.$line.bed|sort -nk 7,7 > $g.$line.bed.tmp && mv -f $g.$line.bed.tmp $g.$line.bed
done < $g.candidateScafs.txt
fi
done
