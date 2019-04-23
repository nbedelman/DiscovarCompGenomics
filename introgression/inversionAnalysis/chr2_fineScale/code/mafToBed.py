#!/usr/bin/env python

#creates bed files for candidate inversion breakpoint scaffolds
#modified from SWD
#updated Jan 8, 2018
#Nate Edelman

import sys
from operator import itemgetter

mafAlignment=sys.argv[1]
outFile=sys.argv[2]

class mafAttributes(object):
	def __init__(self,mafLine):
		line=mafLine.split()
		try:
			self.type=line[0]
			self.label=line[1]
			try:
				self.start=int(line[2])
			except ValueError:
				self.start=0
			try:
				self.end=self.start+int(line[3])
			except ValueError:
				self.end=0
			try:
				self.strand=line[4]
			except IndexError:
				self.strand="NA"
			try:
				self.scafLength=int(line[5])
			except IndexError:
				self.scafLength=0
		except IndexError:
			self.type="NA"
			self.label="NA"
			self.start="NA"
			self.end="NA"
			self.length="NA"
			self.strand="NA"
			self.scafLength="NA"

	def getType(self):
		return self.type
	def getLabel(self):
		return self.label
	def getStart(self):
		if self.getStrand()=="+":
			return self.start
		else:
			return self.getScafLength()-self.end+1
	def getEnd(self):
		if self.getStrand()=="+":
			return self.end
		else:
			return self.getScafLength()-self.start+1
	def getLength(self):
		return self.end-self.start
	def getStrand(self):
		return self.strand
	def getScafLength(self):
		return self.scafLength
	def getBeforeOverhang(self):
		return self.getStart()
	def getAfterOverhang(self):
		return self.getScafLength()-self.getEnd()


def mafToBedList(mafFile):
	'''  Takes a maf file (like the output of hal2maf)
	outputs a list. The list contains ordered attributes to write to a bed file'''
	alignmentList=[]
	id=0
	candidateName=''
	blockOpen=False
	plusColor="0,0,255"
	minusColor="245,145,50"
	with open(mafFile, "r") as file:
		#Skip all commented out lines -  there can be many at the beginning
		for line in file:
			if not "#" in line:
				maf=mafAttributes(line)
				if maf.getType()=='s':
					#Use the very first s line to define the candidate name
					if candidateName == '':
						candidateName = maf.getLabel().split(".")[1]
						#For every other line, we can now distinguish between candidate and reference sequences
						#Use the candidate line in a block to define the overhangs, and decide if
					if  candidateName in maf.getLabel():
						blockOpen=True
						length=maf.getLength()
						beforeOverhang=maf.getBeforeOverhang()
						afterOverhang=maf.getAfterOverhang()
						id+=1
						if length < 100:
							blockOpen=False
					#Use the query line to define the query name, position, and strand
					elif blockOpen and (not candidateName in maf.getLabel()):
						scaffold = maf.getLabel().split(".")[1]
						if maf.getStrand()=="+":
							color=plusColor
							thinStart=maf.getStart()-beforeOverhang
							thinEnd=maf.getEnd()+afterOverhang
						else:
							color=minusColor
							thinStart=maf.getStart()-afterOverhang
							thinEnd=maf.getEnd()+beforeOverhang
						thisAlign=[scaffold,thinStart,thinEnd,candidateName+str(id),0,maf.getStrand(),maf.getStart(),maf.getEnd(),color]
						alignmentList.append(thisAlign)
				else:
					blockOpen=False
	return alignmentList



def writeToBed(outFile,bedList):
    f=open(outFile, "w")
    for i in range(len(bedList)):
        toWrite=""
        for j in range(len(bedList[i])-1):
            toWrite+=(str(bedList[i][j])+"\t")
        toWrite+=(str(bedList[i][-1])+"\n")
        f.write(toWrite)
    f.close()


def runAll(mafAlignment,outFile):
	bedList = mafToBedList(mafAlignment)
	writeToBed(outFile,bedList)

runAll(mafAlignment,outFile)
