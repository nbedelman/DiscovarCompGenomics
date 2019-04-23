directory=$1
bedBase=$2
mafDir=$3
fastaDir=$4
targetGenomes=$5
refGenome=$6
halFile=$7

baseNum=0
num=0
mkdir $directory\_sub_$baseNum
thisDir=$directory\_sub_$baseNum
for gene in $(ls "$directory/$bedBase"*)
do cp $gene $directory\_sub_$baseNum
let num="$num+1"
if [ $num == 1000 ];then
sbatch code/extractFastaAligns.slurm $thisDir $mafDir $fastaDir $targetGenomes $refGenome $halFile
let baseNum="$baseNum+1"
thisDir=$directory\_sub_$baseNum
num=0
fi
done
