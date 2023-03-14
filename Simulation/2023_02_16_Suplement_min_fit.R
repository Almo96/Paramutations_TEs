library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)                               
theme_set(theme_bw())

p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

df<-read.table("2023_02_16_Simulation_supp_fit", fill = TRUE, sep = "\t")
names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
             "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
             "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p0_x0", "p0_x0.01", "p10_x0", "p10_x0.01"))

g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase), alpha = 1, linewidth = 0.7)+
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
  facet_wrap(~sampleid, ncol = 2, nrow = 2, labeller = labeller(sampleid = 
                                                                  c("p0_x0" = "para = 0%, x = 0 ",
                                                                    "p0_x0.01" = "para = 0%, x = 0.01",
                                                                    "p10_x0" = "para = 10%, x = 0 ",
                                                                    "p10_x0.01" = "para = 10%, x = 0.01 ")))


plot(g)


df2<-read.table("2023_02_16_Simulation_supp_clu_para", fill = TRUE, sep = "\t")
names(df2)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
             "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
             "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

unique(df2$sampleid)
df2$sampleid <- factor(df2$sampleid, levels=c("c0.05_p0", "c0.05_p10", "c0.5_p0", "c0.5_p10", "c5_p0", "c5_p10"))


df2_minw <- subset(df2, gen == 5000)


g2 <- ggplot(df2_minw, aes(x = sampleid, y = 1-(min_w)))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("c0.05_p0", "c0.05_p10"),
                                 c("c0.5_p0", "c0.5_p10"),
                                 c("c5_p0", "c5_p10")),
              map_signif_level = TRUE)+
  ylab("fitness cost")+
  xlab("")+
  theme(axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=24))


plot(g2)
