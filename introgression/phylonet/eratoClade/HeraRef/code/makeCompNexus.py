#!/usr/bin/env python
import sys
import itertools

netFile=sys.argv[1]
nexFile=sys.argv[2]
calc=sys.argv[3]

def compileNexus(netFile,outFile,calc):
    out=open(outFile, "w")
    out.write('#NEXUS\n\nBEGIN NETWORKS;\n\n')
    nets=open(netFile, "r")
    netIDs=[]
    for line in nets:
        atts=line.split()
        netIDs.append(atts[0])
        out.write('''Network %s = %s \n''' % (atts[0],atts[1]))
    out.write('\n END;\n\nBEGIN PHYLONET;\n\n')
    for comb in itertools.combinations(netIDs,2):
        out.write('''Cmpnets %s %s -m %s;\n''' % (comb[0],comb[1],calc))
    out.write('END;')
    out.close()

compileNexus(netFile,nexFile,calc)
