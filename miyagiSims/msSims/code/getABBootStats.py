import sys
import statistics

inFile=sys.argv[1]
outFile=sys.argv[2]

def getStats(inFile):
    i=open(inFile,"r")
    o=open(outFile,"w")
    fd=[]
    D=[]
    i.readline()
    for line in i:
        atts=line.split(",")
        if atts[1] != "nan":
            D.append(float(atts[1]))
        if float(atts[0])<0:
            pass
        elif atts[0] != "nan":
            fd.append(float(atts[0]))
    print('''fd=%f (%f-%f) ; D=%f (%f-%f)''' % (statistics.mean(fd),max(fd),min(fd),statistics.mean(D),max(D),min(D)))
    o.write('''%f,%f,%f,%f,%f,%f''' % (statistics.mean(fd),max(fd),min(fd),statistics.mean(D),max(D),min(D)))
getStats(inFile)
