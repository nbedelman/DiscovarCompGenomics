#!/usr/bin/env python

#if a run did not converge, take the ending states of all the gene trees and the network, and make a new nexus using them as the starting states.

import sys
import os
from Bio import SeqIO
import subprocess

originalNexus=sys.argv[1]
resultsDir=sys.argv[2]
output=sys.argv[3]

def getPreviousResults(resultsDir):
    geneTreeBlock="BEGIN TREES;\n"
    networkBlock="BEGIN NETWORKS;\n"
    for subdir, dirs, files in os.walk(resultsDir):
        for logFile in files:
            if "tree" in logFile:
                locusName=logFile.split("_")[1].split(".")[0]
                getLastLine='''tail -1 %s''' % (resultsDir+"/"+logFile)
                lastLine=subprocess.check_output(getLastLine, shell=True)
                lastTree=lastLine.split()[1]
                # f=open(resultsDir+"/"+logFile,"r")
                # for line in f:
                #     pass
                # lastTree=line
                geneTreeBlock+='''Tree %s = %s\n''' % (locusName,lastTree.strip())
                # f.close()
            elif "network" in logFile:
                getLastLine='''tail -1 %s''' % (resultsDir+"/"+logFile)
                lastLine=subprocess.check_output(getLastLine, shell=True)
                lastNet=lastLine.split()[1]
                # f=open(resultsDir+"/"+logFile,"r")
                # for line in f:
                #     pass
                # lastTree=line
                networkBlock+='''Network net1 = %s\n''' % (lastNet.strip())
                # f.close()
    geneTreeBlock+="END;\n\n"
    networkBlock+="END;\n\n"
    return [geneTreeBlock,networkBlock]

def parseOriginalNexus(originalNexus):
    startAndLoci=""
    phyloBlock=""
    inPhyloBlock=False
    f=open(originalNexus,"r")
    numLoci=0
    for line in f:
        if not inPhyloBlock:
            if "PHYLONET" in line:
                inPhyloBlock=True
                phyloBlock+=line
            else:
                if 'locus' in line:
                    numLoci+=1
                startAndLoci+=line
        elif "MCMC_SEQ" in line:
            if "sgt" not in line:
                orig=line.split(";")[0:-1]
                locusString=''
                for i in range(numLoci):
                    locusString+='''locus%i,''' % (i)
                new='''%s %s;\n''' % (";".join(orig), "-sgt (all) -snet (net1)")
                new=new.replace("all",locusString[:-1])
                phyloBlock+=new
            else:
                phyloBlock+=line
        else:
            phyloBlock+=line
    f.close()
    return [startAndLoci,phyloBlock]

def writeNewNexus(originalNexus,resultsDir,output):
    previousResults=getPreviousResults(resultsDir)
    geneTreeBlock=previousResults[0]
    networkBlock=previousResults[1]
    originalParts=parseOriginalNexus(originalNexus)
    startAndLoci=originalParts[0]
    phyloBlock=originalParts[1]
    out=open(output,"w")
    out.write(startAndLoci+"\n")
    out.write(geneTreeBlock+"\n")
    out.write(networkBlock+"\n")
    out.write(phyloBlock)
    out.close()

writeNewNexus(originalNexus,resultsDir,output)
