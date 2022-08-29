library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
theme_set(theme_bw())

p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation")

df<-read.table("2022_08_11_Simulation_2_supplementary", fill = TRUE, sep = "\t")
names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
             "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p1_0.1", "p10_0.1", "p100_0.1", "p1_1", "p10_1", "p100_1", "p1_10", "p10_10", "p100_10", "p1_100", "p10_100", "p100_100"))

g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  ylab("insertions per diploid individual")+xlab("generations")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid, ncol = 3, labeller = labeller(sampleid = 
                                                        c("p1_0.1" = "Paramutable loci = 1% & Trigger loci = 0.1%",
                                                          "p10_0.1" = "Paramutable loci = 10% & Trigger loci = 0.1%",
                                                          "p100_0.1" = "Paramutable loci = 100% & Trigger loci = 0.1%",
                                                          "p1_1" = "Paramutable loci = 1% & Trigger loci = 1%",
                                                          "p10_1" = "Paramutable loci = 10% & Trigger loci = 1%",
                                                          "p100_1" = "Paramutable loci = 100% & Trigger loci = 1%",
                                                          "p1_10" = "Paramutable loci = 1% & Trigger loci = 10%",
                                                          "p10_10" = "Paramutable loci = 10% & Trigger loci = 10%",
                                                          "p100_10" = "Paramutable loci = 100% & Trigger loci = 10%",
                                                          "p1_100" = "Paramutable loci = 1% & Trigger loci = 100%",
                                                          "p10_100" = "Paramutable loci = 10% & Trigger loci = 100%",
                                                          "p100_100" = "Paramutable loci = 100% & Trigger loci = 100%")))

plot(g)
