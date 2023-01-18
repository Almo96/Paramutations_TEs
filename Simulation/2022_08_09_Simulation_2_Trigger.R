library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
theme_set(theme_bw())

p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

df<-read.table("2022_08_09_Simulation_2_Trigger", fill = TRUE, sep = "\t")
names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
             "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p30_0.1", "p30_1", "p30_10","p30_100"))


g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generations")+
  ylab("insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  ggtitle("piRNA clusters = 0%, paramutable loci = 30%") +
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                                      c("p30_0.1" = "trigger loci = 0.1%",
                                                        "p30_1" = "Trigger loci = 1%",
                                                        "p30_10" = "Trigger loci = 10%",
                                                        "p30_100" = "Trigger loci = 100%")))

plot(g)



df_2<-read.table("2023_01_17_Simulation_2_Trigger_2", fill = TRUE, sep = "\t")
names(df_2)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
               "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
               "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")


df_2$phase <- factor(df_2$phase, levels=c("rapi", "trig", "shot", "inac"))
df_2$sampleid <- factor(df_2$sampleid, levels=c("p10_0.1", "p10_1", "p10_10","p10_100"))


g_2<-ggplot()+
  geom_line(data=df_2,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generations")+
  ylab("insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  ggtitle("piRNA clusters = 0%, paramutable loci = 10%") +
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("p10_0.1" = "Trigger loci = 0.1%",
                                                "p10_1" = "Trigger loci = 1%",
                                                "p10_10" = "Trigger loci = 10%",
                                                "p10_100" = "Trigger loci = 100%")))

plot(g_2)


df_3<-read.table("2023_01_17_Simulation_2_Trigger_3", fill = TRUE, sep = "\t")
names(df_3)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
               "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
               "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")


df_3$phase <- factor(df_3$phase, levels=c("rapi", "trig", "shot", "inac"))
df_3$sampleid <- factor(df_3$sampleid, levels=c("p3_10_0.1", "p3_10_1", "p3_10_10","p3_10_100"))


g_3<-ggplot()+
  geom_line(data=df_3,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generations")+
  ylab("insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  ggtitle("piRNA clusters = 3%, paramutable loci = 10%") +
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("p3_10_0.1" = "Trigger loci = 0.1%",
                                                "p3_10_1" = "Trigger loci = 1%",
                                                "p3_10_10" = "Trigger loci = 10%",
                                                "p3_10_100" = "Trigger loci = 100%")))

plot(g_3)
