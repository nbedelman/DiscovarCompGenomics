#This is to format a window bed file so that it is compatible with Steven's coordinates (some of the chromosomes are inverted in our version relative to his, and the scaffolds have different names)

windowFile=$1
LDfile=$2
transitionFile=~/Dropbox/Chr2_inversion/data/Herato.transitions.tsv

~/Documents/Mallet_Lab/DiscovarCompGenomics/scripts/eratoScaffoldConversion.py $windowFile $windowFile.vbOrder.tmp $transitionFile
awk 'BEGIN {OFS="\t"};{$2=($2<0)?1:$2}1' $windowFile.vbOrder.tmp | sort -k1,1 -k2,2n > $windowFile.vbOrder.bed
bedtools closest -a $windowFile.vbOrder.bed -b $LDfile -wb|awk 'BEGIN {OFS="\t"};{print $1,$2,$3,$4,$10}'> $windowFile.LD.bed
