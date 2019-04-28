#!/usr/bin/env python

import pylab
from operator import itemgetter

import sys

treeFile=sys.argv[1]
outFile=sys.argv[2]

trees=open(treeFile,"r")



treeList=[]
for line in trees:
    chrom=line.split(",")[0].split("_")[0]
    number=int(line.split(",")[0].split("_")[1])
    treeType=line.split(",")[1].strip()
    treeList.append([chrom,number,treeType])

treeList=sorted(treeList,key=itemgetter(0,1))

treeDict={}
num=1
treeType=''
for i in treeList:
    if i[2]==treeType:
        num+=1
    else:
        try:
            treeDict[treeType].append(num)
        except KeyError:
            treeDict[treeType]=[num,]
        num=1
        treeType=i[2]

o=open(outFile,"w")

for k in treeDict.keys():
    for entry in treeDict[k]:
        o.write('''%s,%s\n''' % (k,entry))
