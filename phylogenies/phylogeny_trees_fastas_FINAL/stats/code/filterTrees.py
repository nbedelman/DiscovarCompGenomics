import sys
from StringIO import StringIO
from Bio import Phylo
from Bio.Phylo.Consensus import _BitString

#for full erato clade

treeFile=sys.argv[1]
outFile=sys.argv[2]


#the input here is a tsv with two columns: region and newick tree.
#The output will be a csv with two columns: region and type.

# function to convert a list of nwk trees into a ETE tree objects
def makeTree(tree):
    tree = tree.rstrip()
    if not tree == 'NA':
    #     tree = Tree(tree)
        tree=Phylo.read(StringIO(tree), 'newick')
        return (tree)
    else:
        return ('NA')

def _bitstrs(tree):
    bitstrs = set()
    term_names = [term.name for term in tree.get_terminals()]
    term_names.sort()
    for clade in tree.get_nonterminals():
        clade_term_names = [term.name for term in clade.get_terminals()]
        boolvals = [name in clade_term_names for name in term_names]
        bitstr = _BitString(''.join(map(str, map(int, boolvals))))
        bitstrs.add(bitstr)
    return bitstrs

def compare(tree1, tree2):
    term_names1 = [term.name for term in tree1.get_terminals()]
    term_names2 = [term.name for term in tree2.get_terminals()]
    # false if terminals are not the same
    if set(term_names1) != set(term_names2):
        return False
    # true if _BitStrings are the same
    if _bitstrs(tree1) == _bitstrs(tree2):
        return True
    else:
        return False


t=open(treeFile, "r")
out=open(outFile, "a+")

tree1 = Phylo.read(StringIO('Etal,(((Hsar,Hdem),(Hhsa,Htel)),(Hhim,HeraRef));'), 'newick')
tree1.name='tree1'
tree2 = Phylo.read(StringIO('Etal,(((Hsar,Hdem),Htel),(Hhsa,(Hhim,HeraRef)));'), 'newick')
tree2.name='tree2'
tree3 = Phylo.read(StringIO('Etal,(((Hsar,Hdem),Hhsa),(Htel,(Hhim,HeraRef)));'), 'newick')
tree3.name='tree3'
tree4 = Phylo.read(StringIO('Etal,((Hsar,Hdem),(Htel,(Hhsa,(Hhim,HeraRef))));'), 'newick')
tree4.name='tree4/12'
tree5 = Phylo.read(StringIO('Etal,((Hsar,(Htel,Hhsa)),(Hdem,(Hhim,HeraRef)));'), 'newick')
tree5.name='tree5/13'
tree6 = Phylo.read(StringIO('Etal,((Hsar,Hhsa),((Hdem,Htel),(Hhim,HeraRef)));'), 'newick')
tree6.name='tree6'
tree7 = Phylo.read(StringIO('Etal,((Hsar,Htel),(Hdem,(Hhsa,(Hhim,HeraRef))));'), 'newick')
tree7.name='tree7'
tree8 = Phylo.read(StringIO('Etal,(Hsar,((Hdem,Htel),(Hhsa,(Hhim,HeraRef))));'), 'newick')
tree8.name='tree8/11'
tree9 = Phylo.read(StringIO('Etal,((Hsar,(Htel,(Hdem,Hhsa))),(Hhim,HeraRef));'), 'newick')
tree9.name='tree9'
tree10 = Phylo.read(StringIO('Etal,((Hsar,(Hdem,Hhsa)),(Htel,(Hhim,HeraRef)));'), 'newick')
tree10.name='tree10'
tree14 = Phylo.read(StringIO('Etal,((Hsar,Hhsa),(Htel,(Hdem,(Hhim,HeraRef))));'), 'newick')
tree14.name='tree14'
tree15 = Phylo.read(StringIO('Etal,((Hsar,Htel),((Hdem,Hhsa),(Hhim,HeraRef)));'), 'newick')
tree15.name='tree15'
tree16 = Phylo.read(StringIO('Etal,(Hsar,(Htel,((Hdem,Hhsa),(Hhim,HeraRef))));'), 'newick')
tree16.name='tree16'
tree17 = Phylo.read(StringIO('Etal,((Hsar,Hdem),((Hhsa,Htel),(Hhim,HeraRef)));'), 'newick')
tree17.name='tree17'

trees=[tree1,tree2,tree3,tree4,tree5,tree6,tree7,tree8,tree9,tree10,tree14,tree15,tree16,tree17]

for entry in t:
    atts=entry.split()
    region=entry.split()[0]
    if len(atts)!=2:
        out.write('''%s,%s\n''' % (region,"NA"))
    else:
        region=entry.split()[0]
        tree=makeTree(atts[1])
        if tree != 'NA':
            for i in tree.get_terminals():
                if i.name=='Etal':
                    Etal=i
            tree.root_with_outgroup(Etal)
            matchName=''
            for i in trees:
                if compare(tree, i):
                    matchName=i.name
                    out.write('''%s,%s\n''' % (region,matchName))
                    break
            if matchName=='':
                out.write('''%s,%s\n''' % (region,"other"))
        else:
        	out.write('''%s,%s\n''' % (region,"NA"))
out.close()
t.close()
