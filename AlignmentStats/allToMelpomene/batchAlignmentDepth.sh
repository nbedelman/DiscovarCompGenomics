scafList=$(cat $1)
halFile=$2
refGenome=$3
output=$4
outputDir=$5


for scaf in $scafList;
do echo $scaf;
sbatch getAlignmentDepth_specSeqs.slurm $halFile $refGenome $output $outputDir $scaf;
sleep 1;
done
