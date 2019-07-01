library(ggplot2)
library(dplyr)

ANCIENT_TEST_1KB="ancient_tests.allSimulations.results"
ANCIENT_TEST_5KB="5KB_ancient_tests.allSimulations.results"
DEMO_TEST="demo_tests.allSimulations.results"
RECOMB_TEST="recomb_tests.allSimulations.results"
MS_TEST="msSims/msSims.allSimulations.results"
HSATEL_TEST="hsaTelSims/preliminaryResults.results"
HSATEL_NOR_TEST="hsaTelSims_noRecomb/hsaTel.allSimulations.results"


out="hsaTel_noRecomb"

dat_all <- read.csv(HSATEL_NOR_TEST)
dat_all$trueLambda <- 2*dat_all$popSize*3e-9
dat_sub <- subset(dat_all, popSize==2e6)

adjFrac <- c()
for (line in seq(1,nrow(dat_all))){
  l=dat_all[line,]
  if (l$test == "src_dec"){
    if (l$popSize == 2000000 ){adjFrac <- c(adjFrac,l$frac*.865)
    } else if (l$popSize == 1000000){adjFrac <- c(adjFrac,l$frac*.981)
    } else {adjFrac <- c(adjFrac, l$frac)
    }
  }
  else if (l$test == "src_inc"){
    if (l$popSize == 2000000 ){adjFrac <- c(adjFrac,l$frac*.393)
    } else if (l$popSize == 1000000){adjFrac <- c(adjFrac,l$frac*.632)
    } else if (l$popSize == 500000){adjFrac <- c(adjFrac,l$frac*.865)
    } else {adjFrac <- c(adjFrac, l$frac*.981)
    }
  }
  else {
  if (l$popSize == 2000000 ){adjFrac <- c(adjFrac,l$frac*.632)
  } else if (l$popSize == 1000000){adjFrac <- c(adjFrac,l$frac*.865)
  } else if (l$popSize == 500000){adjFrac <- c(adjFrac,l$frac*.981)
  } else {adjFrac <- c(adjFrac, l$frac)
  }
  }
  }
dat_all$adjFrac <- adjFrac

##plot all tests
fdPlot <- ggplot(data=dat_all)+
  geom_point(aes(x=frac,y=fd,col=as.factor(popSize)))+
  geom_errorbar(aes(x=frac,ymin=fdMin,ymax=fdMax,col=as.factor(popSize)))+
  geom_line(aes(x=frac,y=fd,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  xlab("Simulated introgression fraction")+
  geom_abline(slope=1,intercept=0, lty=2, alpha=.75)+
  scale_color_brewer(type = "qual",palette = 2)
fdPlot
ggsave(fdPlot,file=paste0("plots/fd_",out,".pdf"),device="pdf",width=10, height=7)

fdPlot_adjFrac <- ggplot(data=dat_all)+
  geom_point(aes(x=adjFrac,y=fd,col=as.factor(popSize)))+
  geom_errorbar(aes(x=adjFrac,ymin=fdMin,ymax=fdMax,col=as.factor(popSize)))+
  geom_line(aes(x=adjFrac,y=fd,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  xlab("Simulated introgression fraction")+
  geom_abline(slope=1,intercept=0, lty=2, alpha=.75)+
  scale_color_brewer(type = "qual",palette = 2)
fdPlot_adjFrac
ggsave(fdPlot_adjFrac,file=paste0("plots/fd_ajd_",out,".pdf"),device="pdf",width=10, height=7)

DPlot <- ggplot(data=dat_all)+
  geom_point(aes(x=frac,y=D,col=as.factor(popSize)))+
  geom_errorbar(aes(x=frac,ymin=DMin,ymax=DMax,col=as.factor(popSize)))+
  geom_line(aes(x=frac,y=D,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  xlab("Simulated introgression fraction")+
  geom_abline(slope=1,intercept=0, lty=2, alpha=.75)+
  scale_color_brewer(type = "qual",palette = 2)
DPlot
ggsave(DPlot,file=paste0("plots/D_",out,".pdf"),device="pdf",width=10, height=7)

DPlot_adjFrac <- ggplot(data=dat_all)+
  geom_point(aes(x=adjFrac,y=D,col=as.factor(popSize)))+
  geom_errorbar(aes(x=adjFrac,ymin=DMin,ymax=DMax,col=as.factor(popSize)))+
  geom_line(aes(x=adjFrac,y=D,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  xlab("Simulated introgression fraction")+
  geom_abline(slope=1,intercept=0, lty=2, alpha=.75)+
  scale_color_brewer(type = "qual",palette = 2)
DPlot_adjFrac
ggsave(DPlot_adjFrac,file=paste0("plots/D_adj_",out,".pdf"),device="pdf",width=10, height=7)

quiblIntroPlot <- ggplot(data=dat_all)+
  geom_point(aes(x=frac,y =introProp,col=as.factor(popSize)))+
  geom_line(aes(x=frac,y=introProp,col=as.factor(popSize)))+
  geom_errorbar(aes(x=frac,ymin=minIntroProp,ymax=maxIntroProp,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  xlab("Simulated introgression fraction")+
  geom_abline(slope=1,intercept=0, lty=2, alpha=.75)+
  scale_color_brewer(type = "qual",palette = 2)
quiblIntroPlot
ggsave(quiblIntroPlot,file=paste0("plots/introProp_",out,".pdf"),device="pdf",width=10, height=7)

quiblIntroPlot_adjFrac <- ggplot(data=dat_all)+
  geom_point(aes(x=adjFrac,y =introProp,col=as.factor(popSize)))+
  geom_line(aes(x=adjFrac,y=introProp,col=as.factor(popSize)))+
  geom_errorbar(aes(x=adjFrac,ymin=minIntroProp,ymax=maxIntroProp,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  ylim(c(0,1))+
  xlab("Simulated introgression fraction")+
  geom_abline(slope=1,intercept=0, lty=2, alpha=.75)+
  scale_color_brewer(type = "qual",palette = 2)
quiblIntroPlot_adjFrac
ggsave(quiblIntroPlot_adjFrac,file=paste0("plots/introProp_adj_",out,".pdf"),device="pdf",width=10, height=7)

quiblCPlot <- ggplot(data=dat_all)+
  geom_point(aes(x=frac,y=C,col=as.factor(popSize)))+
  geom_line(aes(x=frac,y=C,col=as.factor(popSize)))+
  geom_errorbar(aes(x=frac,ymin=minC,ymax=maxC,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  #ylim(c(0,1))+
  scale_color_brewer(type = "qual",palette = 2)
quiblCPlot
ggsave(quiblCPlot,file=paste0("plots/C_",out,".pdf"),device="pdf",width=10, height=7)

quiblLambdaPlot <- ggplot(data=dat_all)+
  geom_point(aes(x=frac,y=lambda,col=as.factor(popSize)))+
  geom_line(aes(x=frac,y=lambda,col=as.factor(popSize)))+
  geom_errorbar(aes(x=frac,ymin=minLambda,ymax=maxLambda,col=as.factor(popSize)))+
  facet_wrap(facets=~test)+
  xlim(c(0,1))+
  #ylim(c(0,1))+
  geom_line(aes(x=frac, y=trueLambda, col=as.factor(popSize)),lty=2)+
  scale_color_brewer(type = "qual",palette = 2)
quiblLambdaPlot
ggsave(quiblLambdaPlot,file=paste0("plots/lambda_",out,".pdf"),device="pdf",width=10, height=7)


