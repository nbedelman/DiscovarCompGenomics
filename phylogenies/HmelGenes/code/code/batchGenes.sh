bedFile=$1
outBase=$2

num=0
baseNum=0

while read line;
do gene=$(echo $line|awk '{print $4}')
echo $line> $outBase\_$baseNum\_$num\_$gene.bed
num=$((num + 1));
if [ "$num" = 10000 ] ; then
  echo "10000";
  num=0;
  baseNum=$((baseNum + 1));
fi
done < $bedFile
