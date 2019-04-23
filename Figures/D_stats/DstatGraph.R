library(ggplot2)

Dstats=read.csv("~/Dropbox/heliconius/MANUSCRIPT/Sections/supplementaryABBABABA/allHeliconius.singleCopy.D.csv", header=F,
                col.names=c("out","P3","P2","P1","D","SE","Z","best","BABAs","ABBAs","SNPs"))
#Dstats$realD=(Dstats$ABBAs-Dstats$BABAs)/(Dstats$ABBAs+Dstats$BABAs)
Dstats$pVal=2*pnorm(-abs(Dstats$Z))
bestStats <- subset(Dstats,best=='best' & P1 != 'HeraHhimHyb'& P2 != 'HeraHhimHyb'& P3 != 'HeraHhimHyb' & out == 'Etal')
bestStats2 <- subset(Dstats,best=='best')
#par(mfrow=c(2,2))
ZscoreGraph <- ggplot(data=bestStats)+
  geom_histogram(aes(x=Z), bins=200)+
  labs(x="Z score")+
  theme(text = element_text(size=20))+
  geom_vline(xintercept=c(-3.74835,3.74835), lty=2)
  #xlim(c(-40,40))
  # geom_vline(xintercept=c(-10.502,10.502), lty=2)+
  # geom_vline(xintercept=c(-12.80491,12.80491), lty=2)
ZscoreGraph

DscoreGraph <- ggplot(data=bestStats)+
  geom_histogram(aes(x=D), bins=200)+
  labs(x="D statistic")+
  theme(text = element_text(size=20))+
  xlim(c(-.25,.25))
# geom_vline(xintercept=c(-3,3), col="blue", lwd=2)+
# geom_vline(xintercept=c(-6,6), col="red", lwd=2)
DscoreGraph

PvalGraph <- ggplot(data=bestStats)+
  geom_histogram(aes(x=-log(pVal)), bins=200)+
  labs(x="-log p-values", title="p-value")+
  geom_vline(xintercept=c(-log(8.928571e-05)))
# geom_vline(xintercept=c(-6,6), col="red", lwd=2)
PvalGraph

length(which(bestStats$pVal<.05/364))
length(which(bestStats$pVal<.05/1092))
length(which(abs(bestStats$Z)> 3.81))

eratoDstats=read.csv("introgression/ADMIXTOOLS/allHeliconius_eratoCladeD.best.csv", header=T)

length(which(abs(eratoDstats$Z)> 8.892))
bigD=subset(eratoDstats,abs(Z)>8.892)
