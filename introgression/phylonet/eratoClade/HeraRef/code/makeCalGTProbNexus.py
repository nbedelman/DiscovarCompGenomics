#!/bin/env python

import sys
import os
import random
from Bio import SeqIO

output=sys.argv[1]
network=sys.argv[2]
geneTrees=sys.argv[3]
try:
    additionalOptions=sys.argv[4]
except IndexError:
    additionalOptions=''

def writeNexus(output,network, geneTreeFile, additionalOptions):
    o=open(output,"w")
    o.write('#NEXUS\n\nBEGIN NETWORKS;\n\n')
    o.write('''Network net = %s\n\nEND;\n\n''' % (network))
    o.write('BEGIN TREES;\n\n')
    gt=open(geneTreeFile, "r")
    for num,tree in enumerate(gt):
        o.write('''Tree geneTree%i = %s\n'''%(num,tree.strip()))
    o.write('''\nEND;\n\nBEGIN PHYLONET;\n\nCalGTProb net (all) -pl 32 -o %s;\n\nEND;''' % (additionalOptions))
    o.close()


writeNexus(output,network, geneTrees, additionalOptions)
