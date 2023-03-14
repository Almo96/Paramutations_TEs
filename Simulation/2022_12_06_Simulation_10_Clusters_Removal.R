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
  theme(plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24),
        legend.position="none")+
  scale_colour_manual(values=p)+
  scale_x_continuous(breaks = seq(0, 5000, by = 2500))+
  facet_wrap(~sampleid, ncol = 5, nrow = 2, labeller = labeller(sampleid = 
                                                                  c("p0_0" = "para = 0%, clu -0 ",
                                                                    "p0_1" = "para = 0%, clu -1 ",
                                                                    "p0_2" = "para = 0%, clu -2 ",
                                                                    "p0_3" = "para = 0%, clu -3 ",
                                                                    "p0_4" = "para = 0%, clu -4 ",
                                                                    "p10_0" = "para = 10%, clu -0 ",
                                                                    "p10_1" = "para = 10%, clu -1 ",
                                                                    "p10_2" = "para = 10%, clu -2 ",
                                                                    "p10_3" = "para = 10%, clu -3 ",
                                                                    "p10_4" = "para = 10%, clu -4 ")))


plot(g)


df2<-read.table("2022_12_10_Simulation_10_Clusters_removal_3000", fill = TRUE, sep = "\t")
names(df2)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
              "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
              "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df2$phase <- factor(df2$phase, levels=c("rapi", "trig", "shot", "inac"))
df2$sampleid <- factor(df2$sampleid, levels=c("p0_0", "p0_1", "p0_2", "p0_3", "p0_4", "p10_0", "p10_1", "p10_2", "p10_3", "p10_4"))

g2<-ggplot()+
  geom_line(data=df2,aes(x=gen,y=avtes,group=rep,color=phase), alpha = 1, linewidth = 0.7)+
  geom_vline(xintercept = 3000, linetype="dashed", color = "black", linewidth = 0.7)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  theme(plot.title = element_text(size=24),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24),
        legend.position="none")+
  scale_colour_manual(values=p)+
  scale_x_continuous(breaks = seq(0, 5000, by = 1000))+
  facet_wrap(~sampleid, ncol = 5, nrow = 2, labeller = labeller(sampleid = 
                                                                  c("p0_0" = "para = 0, clu -0 ",
                                                                    "p0_1" = "para = 0, clu -1 ",
                                                                    "p0_2" = "para = 0, clu -2 ",
                                                                    "p0_3" = "para = 0, clu -3 ",
                                                                    "p0_4" = "para = 0, clu -4 ",
                                                                    "p10_0" = "para = 10, clu -0 ",
                                                                    "p10_1" = "para = 10, clu -1 ",
                                                                    "p10_2" = "para = 10, clu -2 ",
                                                                    "p10_3" = "para = 10, clu -3 ",
                                                                    "p10_4" = "para = 10, clu -4 ")))


plot(g2)


df3<-read.table("2022_12_10_Simulation_10_Clusters_removal_4000", fill = TRUE, sep = "\t")
names(df3)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
              "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
              "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df3$phase <- factor(df3$phase, levels=c("rapi", "trig", "shot", "inac"))
df3$sampleid <- factor(df3$sampleid, levels=c("p0_0", "p0_1", "p0_2", "p0_3", "p0_4", "p10_0", "p10_1", "p10_2", "p10_3", "p10_4"))


g3<-ggplot()+
  geom_line(data=df3,aes(x=gen,y=avtes,group=rep,color=phase), alpha = 1, linewidth = 0.7)+
  geom_vline(xintercept = 4000, linetype="dashed", color = "black", linewidth = 0.7)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  scale_x_continuous(breaks = seq(0, 5000, by = 2500))+
  facet_wrap(~sampleid, ncol = 5, nrow = 2, labeller = labeller(sampleid = 
                                                                  c("p0_0" = "para = 0, clu -0 ",
                                                                    "p0_1" = "para = 0, clu -1 ",
                                                                    "p0_2" = "para = 0, clu -2 ",
                                                                    "p0_3" = "para = 0, clu -3 ",
                                                                    "p0_4" = "para = 0, clu -4 ",
                                                                    "p10_0" = "para = 10, clu -0 ",
                                                                    "p10_1" = "para = 10, clu -1 ",
                                                                    "p10_2" = "para = 10, clu -2 ",
                                                                    "p10_3" = "para = 10, clu -3 ",
                                                                    "p10_4" = "para = 10, clu -4 ")))


plot(g3)


dfro<-subset(df2, gen == 5000)
s1<-subset(dfro, sampleid == 'p0_0')
s2<-subset(dfro, sampleid == 'p0_4')
wilcox.test(s1$avtes,s2$avtes)
#Wilcoxon rank sum test with continuity correction
#data:  s1$avtes and s2$avtes
#W = 1388, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0



dfro<-subset(df2, gen == 5000)
s1<-subset(dfro, sampleid == 'p10_0')
s2<-subset(dfro, sampleid == 'p10_4')
wilcox.test(s1$avtes,s2$avtes)
#Wilcoxon rank sum test with continuity correction
#data:  s1$avtes and s2$avtes
#W = 4739, p-value = 0.5244
#alternative hypothesis: true location shift is not equal to 0

