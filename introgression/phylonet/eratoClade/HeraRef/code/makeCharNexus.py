#!/usr/bin/env python
import sys
import os

networks=sys.argv[1]
nexBase=sys.argv[2]
calc=sys.argv[3]

def compileNexus(network,outFile,calc):
    out=open(outFile, "w")
    out.write('''#NEXUS\n\nBEGIN NETWORKS;\n\n\
Network net = %s\n\
END;\n\n\
BEGIN PHYLONET;\n\n\
Charnet net -m %s;\n\n\
END;''' % (network,calc))
    out.close


nets=open(networks,"r")
for num,net in enumerate(nets):
    outFile='''%s_%s_%i.nexus''' % (nexBase,os.path.basename(networks),num)
    compileNexus(net,outFile,calc)
