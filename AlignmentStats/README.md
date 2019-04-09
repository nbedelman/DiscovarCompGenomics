#Alignment Stats  

Here, I just wand to get some basic statistics about what was aligned. I'll compare the alignment of melpomene discovar erato reference to melpomene reference eueides, and bombyx to mel ref. I'll also look at the regions where all 25 of the genomes align.

##Strategy:

for the pairwise comparisons, do halAlignmentDepth and wigToBed to get regions of Alignment

Then, extract exons, introns, and intergenic regions from Hmel2.gff.

To get the right bed files, first extract the subset of the gff file that is actually in the alignment (greater than or equal to 1kb scaffolds). do:

```shell
grep <scaffold identifier> <1kb+ genome>  | awk '{print $1}'|sed 's/>//g' > GE1KB_scaffs.txt
```

then:
```python
scafs=[]
s=open("GE1KB_scaffs.txt","r")
for line in s:
  scafs.append(line.strip())
f=open("gffFile.gff","r")
o=open("newGffFile.gff","w")
for line in f:
  scaf=line.split()[0]
  if scaf in scafs:
    o.write(line)
s.close()
f.close()
o.close()
```


exons:
```shell
awk '$3=="exon" {print}' gffFile.gff > exon.gff
```

introns:
```shell
awk '$3=="gene" {print}' gffFile.gff > gene.gff
subtractBed -a gene.gff -b exon.gff > intron.gff
```

intergenic:
```python
#make full scaffold bed file
from Bio import SeqIO
genome=SeqIO.parse(<genomeFile>,"fasta")
for record in genome:
   entry=[record.id,1,len(record)]
   intervals.append(entry)

out=open(<fullScaffolds.bed>,"w")
for e in intervals:
  out.write(e[0]+"\t"+str(e[1])+"\t"+str(e[2])+"\n")

```
OR, if the gff has region entries:
```shell
awk '$3=="region" {print}' gffFile.gff > fullScaffolds.gff
```


```shell
subtractBed -a fullScaffolds.bed -b gene.gff > intergenic.bed
```

Then, for all types, get alignment with:

```shell
intersectBed -wo -a depthBedFile -b gffFile > outputFile
```
The wo option outputs the number of base pairs involved in the overlap. I'll then compare the number of bp overlapped to the total number of bp in the gff file.
