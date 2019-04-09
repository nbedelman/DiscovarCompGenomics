#this is the master script for getting the alignment stats from a wig file.
#It includes basic stats like number of bases with each alignment depth
#and slightly more complex things like sizes of blocks of certain depth.
#I also want to distinguish between genic and non-genic regions.

export PATH=$PATH:~/Documents/Mallet_Lab/DiscovarCompGenomics/scripts

wigFile=$1
baseFile=$(basename $wigFile .wig)

#first, see what the highest depth is
highDepth=$(grep -v fixed $wigFile | awk '{for(i=1;i<=NF;i++) if($i>maxval) maxval=$i;}; END { print maxval;}' -)

#next, get the most basic statistics on number of bases with certain depths
wigAlignmentStats.py $wigFile
