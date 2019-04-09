def makeTransferDict(NWtoDFFile, DFtoScafFile):
    '''takes the scaffold name transfer file downloaded from ftp://ftp.ncbi.nih.gov/genomes/Bombyx_mori/scaffold_names
    and creates a dictionary with DF names as keys and NW names as values. 
    Then takes scaffold file downloaded from http://sgp.dna.affrc.go.jp/KAIKObase/keyword_search.php?keyword=&field=all&chr_or_scaf=&chr_id=all&chr_start=&chr_end=&scaf_id=&scaf_start=&scaf_end=&graph_view=on&page=16
    and creates a new dictionary with scaf as keys and NW as values'''
    f1=open(NWtoDFFile, "r")
    DF_NWdict={}
    for line in f1:
        atts=line.split()
        DF=atts[3].split(".")[0]
        NW=atts[2].split(".")[0]
        DF_NWdict[DF]=NW
    f1.close()
    
    scaf_NWdict={}
    f2=open(DFtoScafFile, "r")
    for line in f2:
        atts=line.split()
        scaf=atts[0]
        DF=atts[2]
        scaf_NWdict[scaf]=DF_NWdict[DF]
    f2.close()
    return scaf_NWdict
        
def liftoverGFF(gffFile,transferDict):
    genes=open(gffFile,"r")
    output=open(gffFile+"_transferred.gff","a+")
    errs=open("errors.gff","a+")
    for line in genes:
        #print line
        atts=line.split()
        try:
            atts[0]=transferDict[atts[0]]
            writable=''
            for item in atts[:-1]:
                writable+=item+"\t"
            writable+=atts[-1]+"\n"
            output.write(writable)
        except KeyError:
            pass
            errs.write(line)
            output.write(line)
    output.close()
    