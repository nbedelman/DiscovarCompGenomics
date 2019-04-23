#This script produces the densitree in Figure 1A. The colors and background tree were added with Illustrator.
#The data was generated using code found in code/treesInWindows

library(devtools)
source("https://bioconductor.org/biocLite.R")
library(tidyr)
library(dplyr)
library(phytools)
library(phangorn)
library(ggtree)
library(ape)
library(bindrcpp)

setwd(dir = "/Users/nbedelman/Dropbox/HeliconiusComparativeGenomics/")


############## Define Dataset ###############
treeList <- "datasets/heliconiini.10KBAbutting.trees.tsv"
STATS="phylogenies_stats/heliconiini.10KBAbutting.windowStats.tsv"

rootTaxon=c("Etal")
helTaxa <- c("Hcyd","Htim","HmelDisco","HmelRef","Hpar","Hbes","Hnum","Ldor",
             "HeraDisco","HeraRef","HeraHhimHyb","Hhim", "Hhsa","Htel","Hdem","Hsar")
############## define base Tree #################
#This allows us to order the taxa in a reasonable way
baseTree=read.tree(text="(Etal,(((Hsar,Hdem),(Htel,(Hhsa,(Hhim,HeraHhimHyb,(HeraRef,HeraDisco))))),(Ldor,((Hnum,Hbes,Hpar),((HmelRef,HmelDisco),(Htim,Hcyd))))));")
noOut <- read.tree(text="(((Hsar,Hdem),(Htel,(Hhsa,(Hhim,HeraHhimHyb,(HeraRef,HeraDisco))))),(Ldor,((Hnum,Hbes,Hpar),((HmelRef,HmelDisco),(Htim,Hcyd)))));")
############## format and filter data ###############

treeTips=baseTree$tip.label

#stats
stats <- read.delim(STATS)
stats$block_id <- as.character(stats$block_id)
stats$mean_bootstrap_support=as.numeric(as.character(stats$mean_bootstrap_support))

trees <- read.delim(sep=" ",file=treeList, header=F, col.names = c("block_id","tree"))

treeStats <- inner_join(stats,trees)
treeStats$tree <- as.character(treeStats$tree)
treeStats$phylo <- read.tree(text=treeStats$tree)

highQualNewick <- highQualTrees <- subset(treeStats, mean_bootstrap_support >=80 & fullLength >= 2000 & p_missing < .2)$tree
highQualTrees <- subset(treeStats, mean_bootstrap_support >=80 & fullLength >= 2000 & p_missing < .2)$phylo

treeInputs=as.multiPhylo(highQualTrees)


#root the trees with eueides, then remove eueides for visualization
rootedTrees <- root(treeInputs,rootTaxon)

noOutTrees<-lapply(rootedTrees,function(t,tips)
  drop.tip(t,setdiff(t$tip.label,tips)),tips=noOut$tip.label)
class(noOutTrees)<-"multiPhylo"
allTrees <- noOutTrees

#sample 500 trees so they are visible
helSample <- allTrees[sample(1:length(allTrees), 500, replace=FALSE)]
class(helSample)<-"multiPhylo"

#plot trees with normalized scale
densiTree(helSample, type="cladogram", alpha=0.1, consensus=noOut,scaleX=F,direction = "rightwards", col="#4D4D4D")


