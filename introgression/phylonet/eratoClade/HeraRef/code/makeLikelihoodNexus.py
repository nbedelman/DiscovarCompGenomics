#!/bin/env python

import sys
import os
import random
from Bio import SeqIO

nexFile=sys.argv[1]
treeFile=sys.argv[2]
numRetics=sys.argv[3]
outFile=sys.argv[4]
try:
    additionalOptions=sys.argv[5]
except IndexError:
    additionalOptions=''

def writeNexus(nexFile,treeFile, numRetics,outFile,additionalOptions):
    o=open(nexFile,"w")
    o.write('#NEXUS\n\n')
    o.write('BEGIN TREES;\n\n')
    gt=open(treeFile, "r")
    for num,tree in enumerate(gt):
        o.write('''Tree geneTree%i = %s\n'''%(num,tree.split()[1].strip()))
    o.write('''\nEND;\n\nBEGIN PHYLONET;\n\nInferNetwork_ML (all) %s -pl 32 -bl -di %s %s;\n\nEND;''' % (numRetics,additionalOptions,outFile))
    o.close()


writeNexus(nexFile,treeFile, numRetics,outFile,additionalOptions)
