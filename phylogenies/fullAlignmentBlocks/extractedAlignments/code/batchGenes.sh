bedFile=$1
outBase=$2

num=0

while read line;
do
echo $line>$outBase\_$num.bed
num=$((num + 1));
done < $bedFile
