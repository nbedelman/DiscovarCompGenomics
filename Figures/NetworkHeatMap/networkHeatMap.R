

corr=read.delim("~/Documents/Mallet_Lab/DiscovarCompGenomics/introgression/phyloNet/networkComparison/eratoClade.fullAlign.100loc.long.onlyBest.highConf.luay.tsv",
                h=F,col.names=c("net1","net2","distance"), as.is=TRUE)
corr.names <- sort(unique(c(corr$net1, corr$net2)))
# create a matrix of the right size and put names on it
corr.dist <- matrix(0, length(corr.names), length(corr.names))
dimnames(corr.dist) <- list(corr.names, corr.names)
# create indices by converting names to numbers and create the normal and reversed
# to fill in all the matrix
corr.ind <- rbind(cbind(match(corr$net1, corr.names), match(corr$net2, corr.names)),
                    cbind(match(corr$net2, corr.names),match(corr$net1, corr.names)))
corr.dist[corr.ind] <- rep(corr$distance, 2)
#corr.dist
new.palette=colorRampPalette(c("black","red","yellow","white"),space="rgb") 
heatmap(corr.dist,symm=TRUE, col=new.palette(20))


corr.uniq=unique(corr.dist)
transposed <- t(corr.uniq)
trans.uniq <- unique(transposed)
heatmap(trans.uniq,symm=TRUE,col=new.palette(20))
