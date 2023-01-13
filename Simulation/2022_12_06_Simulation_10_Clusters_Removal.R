library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
theme_set(theme_bw())

p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

df<-read.table("2022_12_10_Simulation_10_Clusters_removal", fill = TRUE, sep = "\t")
names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
             "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
             "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p0_0", "p0_1", "p0_2", "p0_3", "p0_4", "p10_0", "p10_1", "p10_2", "p10_3", "p10_4"))


g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase), alpha = 1, linewidth = 0.7)+
  geom_vline(xintercept = 2000, linetype="dashed", color = "black", linewidth = 0.7)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  scale_x_continuous(breaks = seq(0, 5000, by = 2500))+
  facet_wrap(~sampleid, ncol = 5, nrow = 2, labeller = labeller(sampleid = 
                                                                  c("p0_0" = "para = 0, clusters -0 ",
                                                                    "p0_1" = "para = 0, clusters -1 ",
                                                                    "p0_2" = "para = 0, clusters -2 ",
                                                                    "p0_3" = "para = 0, clusters -3 ",
                                                                    "p0_4" = "para = 0, clusters -4 ",
                                                                    "p10_0" = "para = 10, clusters -0 ",
                                                                    "p10_1" = "para = 10, clusters -1 ",
                                                                    "p10_2" = "para = 10, clusters -2 ",
                                                                    "p10_3" = "para = 10, clusters -3 ",
                                                                    "p10_4" = "para = 10, clusters -4 ")))


plot(g)
