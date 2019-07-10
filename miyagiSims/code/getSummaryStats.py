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
        if atts[8] != "nan":
            D.append(float(atts[8]))
        if float(atts[9])<0:
            fd.append(0)
        elif atts[9] != "na":
            fd.append(float(atts[9]))
    print('''fd=%f (%f) ; D=%f (%f)''' % (statistics.mean(fd),statistics.stdev(fd),statistics.mean(D),statistics.stdev(D)))
    o.write('''%f,%f,%f,%f''' % (statistics.mean(fd),statistics.stdev(fd),statistics.mean(D),statistics.stdev(D)))
getStats(inFile)
