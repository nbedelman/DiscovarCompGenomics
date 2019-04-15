#!/usr/bin/env python

#Takes a VCF produced by snp-sites from fastas including all erato clade individuals
#specific to erato clade

import sys
import scipy
from scipy import spatial
import os
import csv

vcfFolder=sys.argv[1]
outFile=sys.argv[2]

def getDistanceMatrix(vcfFile):
    vcf=open(vcfFile, "r")
    matrix=[[],[],[],[],[],[],[]]
    snps=0
    for line in vcf:
        if "length" in line:
            length=int(line.split("=")[-1].split(">")[0])
        if (not "#" in line) and length>2000:
            snps+=1
            atts=line.split()
            matrix[0].append(int(atts[9]))
            matrix[1].append(int(atts[10]))
            matrix[2].append(int(atts[11]))
            matrix[3].append(int(atts[12]))
            matrix[4].append(int(atts[13]))
            matrix[5].append(int(atts[14]))
            matrix[6].append(int(atts[15]))
    dist=list(scipy.spatial.distance.pdist(matrix,"hamming"))
    finalDist=[(i*snps)/length for i in dist]
    region=vcfFile.split("/")[-1][:-4]
    finalDist.insert(0,region)
    return finalDist

def write_csv(file, asequence, header=None):
    fp =  open(file, 'a+')
    a = csv.writer(fp, delimiter=',')
    if header:
        a.writerows(header)
    a.writerows(asequence)
    fp.close()

allRegions=[]
for filename in os.listdir(vcfFolder):
    allRegions.append(getDistanceMatrix(os.path.join(vcfFolder, filename)))

write_csv(outFile,allRegions)
