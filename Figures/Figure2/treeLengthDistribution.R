library(ggplot2)
library(RColorBrewer)

treeDist=read.csv("~/Dropbox/introgression/eratoClade/50KBsliding_toErato/eratoClade.50KBSlide.HmelOut.dist", col.names=c("tree","num"))
colors=brewer.pal(8,"Paired")

outBase="~/Dropbox/introgression/eratoClade/50KBsliding_toErato/treeDist."

getBinWidth <- function(vect){
  return(2*(IQR(vect, na.rm=T)/length(vect)^(1/3)))
}


pdf(file = paste0(outBase,"tree2.pdf"), width=2.1676,height=1.3)
g <- ggplot(data=subset(treeDist, tree=="tree2"))+
  geom_histogram(aes(x=num,fill=tree), binwidth=1, col="black", lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[2]))+
  xlim(c(0,max(subset(treeDist, tree=="tree2")$num)+5+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()

pdf(file = paste0(outBase,"tree4.pdf"), width=2.1676,height=1.3)
g <- ggplot(data=subset(treeDist, tree=="tree4"))+
  geom_histogram(aes(x=num,fill=tree), binwidth=1, col="black",lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[6]))+
  xlim(c(0,max(subset(treeDist, tree=="tree4")$num)+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()

pdf(file = paste0(outBase,"tree0.pad.pdf"), width=2.1676,height=1.3)
#padding inversion blocks for display purposes
padFrame <- data.frame(tree=rep("tree0", 9), num=c(40,40,40,18,18,18,43,43,43))
treeDist <- rbind(treeDist,padFrame)
g <- ggplot(data=subset(treeDist, tree=="tree0"))+
  geom_histogram(aes(x=num,fill=tree), binwidth=1, col="black",lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[3]))+
  xlim(c(0,max(subset(treeDist, tree=="tree0")$num)+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()

pdf(file = paste0(outBase,"tree5.pad.pdf"), width=2.1676,height=1.3)
padFrame <- data.frame(tree=rep("tree5", 3), num=c(29,29,29))
treeDist <- rbind(treeDist,padFrame)
g <- ggplot(data=subset(treeDist, tree=="tree5"))+
  geom_histogram(aes(x=num,fill=tree), binwidth=1, col="black",lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[4]))+
  xlim(c(0,max(subset(treeDist, tree=="tree5")$num)+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()

pdf(file = paste0(outBase,"tree6.pdf"), width=2.1676,height=1.3)
g <- ggplot(data=subset(treeDist, tree=="tree6"))+
  geom_histogram(aes(x=num,fill=tree), binwidth=1, col="black",lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[8]))+
  xlim(c(0,max(subset(treeDist, tree=="tree6")$num)+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()

pdf(file = paste0(outBase,"tree9.pdf"), width=2.1676,height=1.3)
g <- ggplot(data=subset(treeDist, tree=="tree9"))+
  geom_histogram(aes(x=num,fill=tree),binwidth=1, col="black",lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[1]))+
  xlim(c(0,max(subset(treeDist, tree=="tree9")$num)+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()

pdf(file = paste0(outBase,"tree3.pdf"), width=2.1676,height=1.3)
g <- ggplot(data=subset(treeDist, tree=="tree3"))+
  geom_histogram(aes(x=num,fill=tree), binwidth=1, col="black",lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[7]))+
  xlim(c(0,max(subset(treeDist, tree=="tree3")$num)+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()

pdf(file = paste0(outBase,"tree7.pdf"), width=2.1676,height=1.3)
g <- ggplot(data=subset(treeDist, tree=="tree7"))+
  geom_histogram(aes(x=num,fill=tree), binwidth=1, col="black",lwd=.1, show.legend=F)+
  scale_fill_manual(values=c(colors[5]))+
  xlim(c(0,max(subset(treeDist, tree=="tree7")$num)+5))+
  theme(text = element_text(size=9))+
  labs(x="",y="")
g
dev.off()