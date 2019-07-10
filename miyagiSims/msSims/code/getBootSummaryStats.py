import sys
import statistics

inFile=sys.argv[1]
outFile=sys.argv[2]

def getStats(inFile):
    i=open(inFile,"r")
    o=open(outFile,"w")
    C=[]
    lamb=[]
    introProp=[]
    i.readline()
    for line in i:
        atts=line.split(",")
        C.append(float(atts[0]))
        lamb.append(float(atts[1]))
        introProp.append(float(atts[2]))
    print('''C=%f (%f-%f) ; lambda=%f (%f-%f); introProp=%f (%f-%f)''' % (statistics.mean(C),max(C),min(C),statistics.mean(lamb),max(lamb),min(lamb),statistics.mean(introProp),max(introProp),min(introProp)))
    o.write('''%f,%f,%f,%f,%f,%f,%f,%f,%f''' % (statistics.mean(C),max(C),min(C),statistics.mean(lamb),max(lamb),min(lamb),statistics.mean(introProp),max(introProp),min(introProp)))
getStats(inFile)
