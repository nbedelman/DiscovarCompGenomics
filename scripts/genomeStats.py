from Bio import SeqIO
from Bio import Seq
import csv
import os
import sys

#Sample	#Seqs	Min	1st Qu.	Median	Mean	3rd Qu.	Max	Total	n50	n90	n95 #seqs>1kb pct>1kb #seqs>5kb pct>5kb
def getGenomeStats(fastaFile):
    lengths=[]
    ge1KB=[]
    ge5KB=[]
    genome=SeqIO.parse(open(fastaFile, "r"), "fasta")
    for record in genome:
        lengths.append(len(record))
        if len(record)>=1000:
            ge1KB.append(len(record))
            if len(record)>=5000:
                ge5KB.append(len(record))
    n50Lens=[]
    for i in lengths:
        for j in range(i):
            n50Lens.append(i)
    n50Lens=sorted(n50Lens)
    n50=getMedian(n50Lens)
    medianSize=getMedian(lengths)
    numSeqs=len(lengths)
    minSize=min(lengths)
    meanSize=float(sum(lengths))/len(lengths)
    maxSize=max(lengths)
    totSize=sum(lengths)
    firstQ=getMedian(lengths[:int(len(lengths)/2)])
    thirdQ=getMedian(lengths[int(len(lengths)/2):])
    n90=n50Lens[int(len(n50Lens)*.9)]
    n95=n50Lens[int(len(n50Lens)*.95)]
    ge1KBseqs=len(ge1KB)
    ge1KBpct=float(sum(ge1KB))/sum(lengths)
    ge5KBseqs=len(ge5KB)
    ge5KBpct=float(sum(ge5KB))/sum(lengths)
    print '''%s,%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t%s/t''' % (fastaFile.split("/")[-1],numSeqs,\
    minSize,firstQ,medianSize,meanSize,thirdQ,maxSize,totSize,n50,n90,n95,ge1KBseqs,ge1KBpct,ge5KBseqs,ge5KBpct)
    return [fastaFile.split("/")[-1],numSeqs,\
    minSize,firstQ,medianSize,meanSize,thirdQ,maxSize,totSize,n50,n90,n95,ge1KBseqs,ge1KBpct,ge5KBseqs,ge5KBpct]

def getMedian(aList):
    if len(aList)%2==1:
        medianSpot=(len(aList)+1)/2
        median=aList[medianSpot]
    elif len(aList)%2 == 0:
        midOne=len(aList)/2
        midTwo=midOne+1
        median=(aList[midOne]+aList[midTwo])/2
    return median

# allGenomeStats=[]
# rootdir="/Users/nbedelman/Documents/Mallet_Lab/18Genomes/Genomic_Analysis/clippedGenomes"
# for subdir, dirs, files in os.walk(rootdir):
#     for file in files:
#         if '.fasta' in file:
#             allGenomeStats.append(getGenomeStats(os.path.join(subdir,file)))



def write_csv(file, asequence, header=None):
    fp =  open(file, 'w')
    a = csv.writer(fp, delimiter=',')
    if header:
        a.writerows(header)
    a.writerows(asequence)
    fp.close()

stats=getGenomeStats(sys.argv[1])
write_csv(sys.argv[2],[stats,])
