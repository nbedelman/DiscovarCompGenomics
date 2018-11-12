#This script produces the densitree in Figure 1A. The colors and background tree were added with Illustrator.
#The data was generated 

library(devtools)
source("https://bioconductor.org/biocLite.R")
library(phytools)
library(phangorn)
library(ggtree)
library(ape)
library(bindrcpp)

setwd(dir = "HeliconiusComparativeGenomics/")


############## Define Dataset ###############
treeList <- "datasets/10KBAbutting.fullAlign.ge500.nwk"


rootTaxon="Etal"

############## define base Tree #################
#This allows us to order the taxa in a reasonable way
baseTree=read.tree(text="(Etal,(((Hsar,Hdem),(Htel,(Hhsa,(Hhim,HeraHhimHyb,(HeraRef,HeraDisco))))),(Ldor,((Hnum,Hbes,Hpar),((HmelRef,HmelDisco),(Htim,Hcyd))))));")

############## format and filter data ###############

treeTips=baseTree$tip.label


treeInputs=read.tree(treeList, keep.multi = T)

#These trees include all lepidoptera; reduce to just Heliconius and Eueides for ease of visualization
pruned<-lapply(treeInputs,function(t,tips)
  drop.tip(t,setdiff(t$tip.label,tips)),tips=treeTips)
class(pruned)<-"multiPhylo"


#midpoint-root the trees with eueides
allTrees=root(pruned,rootTaxon)
allTrees=midpoint(allTrees)

#sample 500 trees so they are visible
helSample <- allTrees[sample(1:length(allTrees), 500, replace=FALSE)]
class(helSample)<-"multiPhylo"

#plot trees with normalized scale
densiTree(helSample, type="cladogram", alpha=0.1, consensus=baseTree,scaleX=T,direction = "rightwards", col="#4D4D4D")
