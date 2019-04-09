#get pairwise alignment statistics for all the genomes in the final hal alignment to a single genome

mkdir -p code
mkdir data

source /n/sw/progressiveCactus-latest/progressiveCactus/environment
module load bedtools

ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal data/fullAlignment.hal
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/heliconius_melpomene_melpomene_hmel2_core_32_85_1.gff data
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/scripts/wigToBed.py code
ln -s /n/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/Hmel2_fullScaffolds.bed data

chmod u+x code/*
export PATH=code/:$PATH

awk '$3 == "exon" {print $1"\t"$4"\t"$5"\t"$9"\t.\t"$7}' data/heliconius_melpomene_melpomene_hmel2_core_32_85_1.gff > data/HmelExon.bed
awk '$3 == "gene" {print $1"\t"$4"\t"$5"\t"$9"\t.\t"$7}' data/heliconius_melpomene_melpomene_hmel2_core_32_85_1.gff > data/Hmelgene.bed

bedtools subtract -a data/Hmelgene.bed -b data/HmelExon.bed > data/HmelIntron.bed

genomes=$(halStats --genomes data/fullAlignment.hal)

for g in $genomes
do
if [[ ! $g == *"Anc"* ]]; then
  echo $g
  sbatch code/pairwiseAlignmentBlocks.slurm data/fullAlignment.hal HmelRef $g\_toHmelRef.wig $g\_toHmelRef.bed $g data/HmelExon.bed data/HmelIntron.bed data/Hmelgene.bed pairwiseAlignmentStats.out
fi
done
