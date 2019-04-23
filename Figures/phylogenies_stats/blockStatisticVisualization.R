#Script to generate Supplementary Figure XX

############## Import Libraries ###############
library(devtools)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plotly)
library(tidyr)
library(knitr)
library(gridExtra)
library(grid)

setwd("/Users/nbedelman/Dropbox/HeliconiusComparativeGenomics/")

########### Define Data Files #############
STATS="phylogenies_stats/10KB_eratoClade/eratoClade_10KBAbutting_HeraRef_windowStats.txt"
STATS2="phylogenies_stats/10KB_eratoClade/heliconiini_HeraRef_fullAlignedSites.windowStats.txt"
out="phylogenies_stats/10KB_eratoClade/eratoClade_10KBAbutting_usedForRecombTests.txt"


####### read in data ######
#stats
stats <- read.delim(STATS)
stats$block_id <- as.character(stats$block_id)
stats$mean_bootstrap_support=as.numeric(as.character(stats$mean_bootstrap_support))

###this is just because we added an additional statistic (fully aligned sites)
stats2 <- read.delim(STATS2)
stats2 <- stats2[,c("block_id","fullAlignedSites")]

stats <- inner_join(stats,stats2)
############### Tidy Data and make appropriate subset #########
stats <- stats %>% drop_na
stats$phi_p <- as.numeric(as.character(stats$phi_p))
stats$transformedP <- -log10(stats$phi_p)
#densitree filter
#stats <- subset(stats, mean_bootstrap_support >=80 & fullLength>=2000 & p_missing<.2)
#eratoClade recomb test filter
stats <- subset(stats,fullLength>=2000 & p_missing<.2)

########### make graphs ########
lengths <- ggplot(data=stats)+
  geom_histogram(aes(x=fullLength))+
  labs(x="Alignment length", y="Number of windows")

pctMissing <- ggplot(data=stats)+
  geom_histogram(aes(x=p_missing))+
  labs(x="Fraction missing", y="Number of windows")

piSites <- ggplot(data=stats)+
  geom_histogram(aes(x=n_pi_sites))+
  labs(x="Phylogenetically informative sites", y="Number of windows")

bootstrap <- ggplot(data=stats)+
  geom_histogram(aes(x=mean_bootstrap_support))+
  labs(x="Mean bootstrap support", y="Number of windows")

breakpoint <- ggplot(data=stats)+
  geom_histogram(aes(x=transformedP))+
  labs(x="Probability of recombination within block", y="Number of windows")

fullAlign <- ggplot(data=stats)+
  geom_histogram(aes(x=fullAlignedSites))+
  labs(x="Fully Aligned Sites", y="Number of windows")

ng <- nullGrob()

grid.arrange(lengths, pctMissing, 
             piSites,bootstrap,
             breakpoint,fullAlign)


g <- arrangeGrob(lengths, pctMissing, 
                 piSites,bootstrap,
                 breakpoint,fullAlign)


ggsave(filename = paste0(out,"individualStats.pdf"),g, width = 10, height=10)
#Pairwise Graphs

oneOne <- ggplot(data=stats)+
  geom_point(aes(x=fullLength,y=n_pi_sites),alpha=.05)+
  labs(y="Phylogenetically informative sites",x="")

oneTwo<- ggplot(data=stats)+
  geom_point(aes(x=fullLength,y=p_missing),alpha=.05)+
  labs(y="Fraction missing",x="")

oneThree<- ggplot(data=stats)+
  geom_point(aes(x=fullLength,y=transformedP),alpha=.05)+
  labs(y="Probability of recombination within block",x="")

oneFour <- ggplot(data=stats)+
  geom_point(aes(x=fullLength,y=mean_bootstrap_support),alpha=.05)+
  labs(y="Mean bootstrap support",x="")

oneFive <- ggplot(data=stats)+
  geom_point(aes(x=fullLength,y=fullAlignedSites),alpha=.05)+
  labs(y="Fully Aligned Sites",x="Alignment length")

twoTwo <- ggplot(data=stats)+
  geom_point(aes(x=n_pi_sites,y=p_missing),alpha=.05)+
  labs(y="",x="")

twoThree <- ggplot(data=stats)+
  geom_point(aes(x=n_pi_sites,y=transformedP),alpha=.05)+
  labs(y="",x="")

twoFour <- ggplot(data=stats)+
  geom_point(aes(x=n_pi_sites,y=mean_bootstrap_support),alpha=.05)+
  labs(y="",x="")

twoFive <- ggplot(data=stats)+
  geom_point(aes(x=n_pi_sites,y=fullAlignedSites),alpha=.05)+
  labs(y="",x="Phylogenetically informative sites")

threeThree <- ggplot(data=stats)+
  geom_point(aes(x=p_missing,y=transformedP),alpha=.05)+
  labs(y="",x="")

threeFour <- ggplot(data=stats)+
  geom_point(aes(x=p_missing,y=mean_bootstrap_support),alpha=.05)+
  labs(y="",x="")

threeFive <- ggplot(data=stats)+
  geom_point(aes(x=p_missing,y=fullAlignedSites),alpha=.05)+
  labs(y="",x="Fraction missing")

fourFour <- ggplot(data=stats)+
  geom_point(aes(x=transformedP,y=mean_bootstrap_support),alpha=.05)+
  labs(y="",x="")

fourFive <- ggplot(data=stats)+
  geom_point(aes(x=transformedP,y=fullAlignedSites),alpha=.05)+
  labs(y="",x="Probability of recombination within block")

fiveFive <- ggplot(data=stats)+
  geom_point(aes(x=mean_bootstrap_support,y=fullAlignedSites),alpha=.05)+
  labs(y="",x="Mean Bootstrap Support")

ng <- nullGrob()
grid.arrange(oneOne, ng, ng,ng,ng,        
             oneTwo, twoTwo,ng,ng,ng,
             oneThree,twoThree,threeThree,ng,ng,
             oneFour,twoFour,threeFour,fourFour,ng,
             oneFive,twoFive,threeFive,fourFive,fiveFive)

g <- arrangeGrob(oneOne, ng, ng,ng,ng,        
                 oneTwo, twoTwo,ng,ng,ng,
                 oneThree,twoThree,threeThree,ng,ng,
                 oneFour,twoFour,threeFour,fourFour,ng,
                 oneFive,twoFive,threeFive,fourFive,fiveFive)

ggsave(filename = paste0(out,"allByAll.jpg"),g, width = 13, height=13)

