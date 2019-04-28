#!/usr/bin/env python

import sys

miyagiOut=sys.argv[1]
converted=sys.argv[2]

def convertSummary(miyagiOut,converted):
    m=open(miyagiOut,"r")
    o=open(converted,"w")
    miyagiString=m.readline()
    allEntries=miyagiString.split("(")
    for entry in allEntries:
        try:
            entry=entry.replace("(","").replace(")","").replace("[","").replace("]","").replace("'","").replace(" ","")
            atts=entry.split(",")
            #print(atts)
            triplet=sorted(atts[0:3])
            outgroup=atts[3]
            C1=float(atts[4])
            C2=float(atts[5])
            prop1=float(atts[6])
            prop2=float(atts[7])
            try:
                numTrees=int(atts[-3])
                BIC=float(atts[-5])
                maxL=float(atts[-6])
            except ValueError:
                numTrees=int(atts[-2])
                BIC=float(atts[-4])
                BIC=float(atts[-5])
            outString='''%s_%s_%s,%s,%f,%f,%f,%f,%i,%f,%f\n''' % (triplet[0],triplet[1],triplet[2],outgroup,C1,C2,prop1,prop2,numTrees,BIC,maxL)
            o.write(outString)
        except IndexError:
            pass
        except ValueError:
            pass

convertSummary(miyagiOut,converted)
