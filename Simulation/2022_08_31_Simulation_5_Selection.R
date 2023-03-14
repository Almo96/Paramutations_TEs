library(tidyverse)
library(ggplot2)

p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")


df<-read.table("2022_08_31_Simulation_5_Selection", fill = TRUE, sep = "\t")
names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed",
             "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
             "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid", "extra")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))

df_x<-subset(df, grepl("^.+(0|1)$", sampleid))
df_noxclu<-subset(df, grepl("^.+(_n)$", sampleid))

g_x<-ggplot()+
  geom_line(data=df_x,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  ylab("TEs insertions per diploid individual")+xlab("generation")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid)

plot(g_x)

g_noxclu<-ggplot()+
  geom_line(data=df_noxclu,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  ylab("TEs insertions per diploid individual")+xlab("generation")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid)

plot(g_noxclu)


g_x_2<-ggplot()+
  geom_line(data=df_x,aes(x=gen,y=(avcli),group=rep,color=phase),alpha=1,size=0.7)+
  ylab("cluster insertions per diploid individual")+xlab("generation")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid)

plot(g_x_2)

g_noxclu<-ggplot()+
  geom_line(data=df_noxclu,aes(x=gen,y=avcli,group=rep,color=phase),alpha=1,size=0.7)+
  ylab("cluster insertions per diploid individual")+xlab("generation")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid)

plot(g_noxclu)


df_x_2 <- subset(df_x, phase %in% c("shot", "inac"))
df_x_3 <- data.frame()
repcheck = 1
x = 1
y = 1

while (x<nrow(df_x_2)) {
  if (repcheck != df_x_2[x, 1]){
    y = 1
  }
  if (y == 1){
    if(df_x_2[x, 12]  == "shot"){
      df_x_3<-rbind(df_x_3,df_x_2[x,])
      y = 2
      repcheck = df_x_2[x, 1]
    }
  }
  if (y == 2){
    if(df_x_2[x, 12] == "inac"){
      df_x_3<-rbind(df_x_3,df_x_2[x,])
      y = 1
    }
  }
  x = x+1
}


e_x <- df_x_3 %>% 
  group_by(sampleid, phase) %>% 
  summarize(mean_avtes = mean(avtes), sd_avtes = sd(avtes), mean_avcli = mean(avcli), sd_avcli = sd(avcli), mean_avpar = mean(avpar), sd_avpar = sd(avpar))

e_x <- e_x[-c(17),]


g_x_3 <- ggplot(e_x, aes(x=phase, y=mean_avtes, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avtes-sd_avtes, ymax=mean_avtes+sd_avtes), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("TEs insertions per diploid individual")+
  xlab("phase")+
  scale_fill_manual(values = c("yellow", "red"))+
  facet_wrap(~sampleid, ncol=4)

plot(g_x_3)

g_x_4 <- ggplot(e_x, aes(x=phase, y=mean_avcli, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avcli-sd_avcli, ymax=mean_avcli+sd_avcli), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("cluster insertions per diploid individual")+
  xlab("phase")+
  scale_fill_manual(values = c("yellow", "red"))+
  facet_wrap(~sampleid, ncol=4)

plot(g_x_4)


df_x_4 <- subset(df_noxclu, phase %in% c("shot", "inac"))
df_x_5 <- data.frame()
repcheck = 1
x = 1
y = 1
while (x<nrow(df_x_4)) {
  if (repcheck != df_x_4[x, 1]){
    y = 1
  }
  if (y == 1){
    if(df_x_4[x, 12]  == "shot"){
      df_x_5<-rbind(df_x_5,df_x_4[x,])
      y = 2
      repcheck = df_x_4[x, 1]
    }
  }
  if (y == 2){
    if(df_x_4[x, 12] == "inac"){
      df_x_5<-rbind(df_x_5,df_x_4[x,])
      y = 1
    }
  }
  x = x+1
}


e_noxclu <- df_x_5 %>% 
  group_by(sampleid, phase) %>% 
  summarize(mean_avcli = mean(avcli), sd_avcli = sd(avcli), mean_avpar = mean(avpar), sd_avpar = sd(avpar))

e_noxclu <- e_noxclu[-c(5,6,11,12,13,14,19,20,25,26),]


g_noxclu_4 <- ggplot(e_noxclu, aes(x=phase, y=mean_avcli, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avcli-sd_avcli, ymax=mean_avcli+sd_avcli), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("cluster insertions per diploid individual")+
  xlab("phase")+
  scale_fill_manual(values = c("yellow", "red"))+
  facet_wrap(~sampleid, ncol=4)

plot(g_noxclu_4)


g_x_5 <- ggplot(e_x, aes(x=phase, y=mean_avpar, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avpar-sd_avpar, ymax=mean_avpar+sd_avpar), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Paramutable loci insertions per diploid individual")+
  xlab("phase")+
  scale_fill_manual(values = c("yellow", "red"))+
  facet_wrap(~sampleid, ncol=4)

plot(g_x_5)

g_noxclu_5 <- ggplot(e_noxclu, aes(x=phase, y=mean_avpar, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avpar-sd_avpar, ymax=mean_avpar+sd_avpar), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Paramutable loci insertions per diploid individual")+
  xlab("phase")+
  scale_fill_manual(values = c("yellow", "red"))+
  facet_wrap(~sampleid, ncol=4)

plot(g_noxclu_5)



df_ori<-read.table("2023_03_05_Simulation_supp_ori", fill = TRUE, sep = "\t")
names(df_ori)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed",
             "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
             "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid", "extra")

df_ori$phase <- factor(df_ori$phase, levels=c("rapi", "trig", "shot", "inac"))
df_ori2 <- data.frame()


#new dataframe with only the first shotgun & the first inactive phase of each replicate
repcheck = 1
x = 1
y = 1
while (x<nrow(df_ori)) {
  if (repcheck != df_ori[x, 1]){
    y = 1
  }
  if (y == 1){
    if(df_ori[x, 12]  == "shot"){
      df_ori2<-rbind(df_ori2,df_ori[x,])
      y = 2
      repcheck = df_ori[x, 1]
    }
  }
  if (y == 2){
    if(df_ori[x, 12] == "inac"){
      df_ori2<-rbind(df_ori2,df_ori[x,])
      y = 1
    }
  }
  x = x+1
}

df_ori2_shot <- subset(df_ori2, phase == "shot")
df_ori2_inac <- subset(df_ori2, phase == "inac")



ori_shot <- ggplot(df_ori2_shot, aes(x = sampleid, y = gen))+
  geom_boxplot()+
  xlab("selection coefficient") +
  ylab("number of origins")

plot(ori_shot)

ori_inac <- ggplot(df_ori2_inac, aes(x = sampleid, y = piori))+
  geom_boxplot()+
  xlab("selection coefficient") +
  ylab("number of origins")

plot(ori_inac)
