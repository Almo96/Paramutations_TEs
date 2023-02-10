library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)
theme_set(theme_bw())

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

p<-c("grey","#1a9850","#ffd700","#d73027")

df<-read.table("2022_08_29_Simulation_3_Clusters", fill = TRUE, sep = "\t")
names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed",
             "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
             "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid", "extra")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p1", "p3", "p10","p50", "p1_10", "p3_10", "p10_10","p50_10"))

df <- df[-c(ncol(df))]

g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  ylab("TEs insertions per diploid individual")+xlab("generation")+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid, ncol=4)

plot(g)


g2<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avcli,group=rep,color=phase),alpha=1,size=0.7)+
  ylab("Cluster insertions per diploid individual")+xlab("generation")+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid, ncol=4)

plot(g2)


df1 <- subset(df, phase %in% c("shot", "inac"))

df2 <- data.frame()

repcheck = 1
x = 1
y = 1
while (x<nrow(df1)+1) {
  if (repcheck != df1[x, 1]){
    y = 1
  }
  if (y == 1){
    if(df1[x, 12]  == "shot"){
      df2<-rbind(df2,df1[x,])
      y = 2
      repcheck = df1[x, 1]
    }
  }
  if (y == 2){
    if(df1[x, 12] == "inac"){
      df2<-rbind(df2,df1[x,])
      y = 1
    }
  }
  x = x+1
}

e <- df2 %>% 
  dplyr::group_by(sampleid, phase) %>% 
  dplyr::summarize(mean_avcli = mean(avcli), sd_avcli = sd(avcli),
            mean_fwpar_yespi = mean(fwpar_yespi))

g2_2 <- ggplot(e, aes(x=phase, y=mean_avcli, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avcli-sd_avcli, ymax=mean_avcli+sd_avcli), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Cluster insertions per diploid individual")+
  xlab("Phase")+
  theme(legend.position="none")+
  scale_y_continuous(limits = c(-0.5, 12.5))+
  scale_fill_manual(values = c("#ffd700", "#d73027"))+
  facet_wrap(~sampleid, ncol=8,labeller = labeller(sampleid = 
                                                    c("p1" = "clu=1%",
                                                      "p3" = "clu=3%",
                                                      "p10" = "clu=10%",
                                                      "p50" = "clu=50%",
                                                      "p1_10" = "clu=1% para=10%",
                                                      "p3_10" = "clu=3% para=10%",
                                                      "p10_10" = "clu=10% para=10%",
                                                      "p50_10" = "clu=50% para=10%")))

plot(g2_2)

dfonlypara <- df %>%
  dplyr::filter(sampleid == "p1_10" | sampleid == "p3_10" | sampleid == "p10_10" | sampleid == "p50_10")


g3<-ggplot()+
  geom_line(data=dfonlypara,aes(x=gen,y=avpar,group=rep,color=phase),alpha=1,size=0.7)+
  ylab("Paramutable site insertions per diploid individual")+xlab("generation")+
  theme(legend.position="none")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_colour_manual(values=p)+
  facet_wrap(~sampleid, ncol=4)

plot(g3)

e_2 <- dfonlypara %>% 
  dplyr::group_by(sampleid, phase) %>% 
  dplyr::summarize(mean_avpar = mean(avpar), sd_avpar = sd(avpar))


e_2 <-subset(e_2, phase!="trig" & phase!="rapi")

g3_2 <- ggplot(e_2, aes(x=phase, y=mean_avpar, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avpar-sd_avpar, ymax=mean_avpar+sd_avpar), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Paramutable loci insertions")+
  xlab("Phase")+
  theme(legend.position="none")+
  scale_y_continuous(limits = c(-0.5, 12.5))+
  scale_fill_manual(values = c("#ffd700", "#d73027"))+
  facet_wrap(~sampleid, ncol=4,labeller = labeller(sampleid =
                                                     c("p1_10" = "clu=1% para=10%",
                                                       "p3_10" = "clu=3% para=10%",
                                                       "p10_10" = "clu=10% para=10%",
                                                       "p50_10" = "clu=50% para=10%")))

plot(g3_2)


e_3 <- df2 %>% 
  dplyr::group_by(sampleid, phase) %>% 
  dplyr::summarize(mean_avpar = mean(avpar), sd_avpar = sd(avpar))


e$clusters<-c("1% piRNA clusters", "1% piRNA clusters", "3% piRNA clusters", "3% piRNA clusters",  
              "10% piRNA clusters", "10% piRNA clusters", "50% piRNA clusters", "50% piRNA clusters",
              "1% piRNA clusters", "1% piRNA clusters", "3% piRNA clusters", "3% piRNA clusters",  
              "10% piRNA clusters", "10% piRNA clusters", "50% piRNA clusters", "50% piRNA clusters")

e$tag<-c("shot", "inac", "shot", "inac", "shot", "inac", "shot", "inac", "shot_para", "inac_para", "shot_para", "inac_para", "shot_para", "inac_para", "shot_para", "inac_para")
e$tag <- factor(e$tag,levels = c("shot", "shot_para", "inac", "inac_para"))
e$clusters <- factor(e$clusters,levels = c("1% piRNA clusters", "3% piRNA clusters", "10% piRNA clusters", "50% piRNA clusters"))

g2_3 <- ggplot(e, aes(x=tag, y=mean_avcli, fill = phase))+ 
  geom_bar(stat = "identity")+
  geom_errorbar( aes(x=tag, ymin=mean_avcli-sd_avcli, ymax=mean_avcli+sd_avcli), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Cluster insertions per diploid individual")+
  xlab("Phase")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_fill_manual(values = c("#ffd700", "#d73027"))+
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~clusters, ncol=4)

plot(g2_3)


e_3$clusters<-c("1% piRNA clusters", "1% piRNA clusters", "3% piRNA clusters", "3% piRNA clusters",  
              "10% piRNA clusters", "10% piRNA clusters", "50% piRNA clusters", "50% piRNA clusters",
              "1% piRNA clusters", "1% piRNA clusters", "3% piRNA clusters", "3% piRNA clusters",  
              "10% piRNA clusters", "10% piRNA clusters", "50% piRNA clusters", "50% piRNA clusters")

e_3$tag<-c("shot", "inac", "shot", "inac", "shot", "inac", "shot", "inac", "shot_para", "inac_para", "shot_para", "inac_para", "shot_para", "inac_para", "shot_para", "inac_para")
e_3$tag <- factor(e_3$tag,levels = c("shot", "shot_para", "inac", "inac_para"))
e_3$clusters <- factor(e_3$clusters,levels = c("1% piRNA clusters", 
"3% piRNA clusters", "10% piRNA clusters", "50% piRNA clusters"))

g3_3 <- ggplot(e_3, aes(x=tag, y=mean_avpar, fill = phase))+ 
  geom_bar(stat = "identity")+
  geom_errorbar( aes(x=tag, ymin=mean_avpar-sd_avpar, ymax=mean_avpar+sd_avpar), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Number of insertions in paramutable loci per diploid individual")+
  xlab("Phase")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  scale_fill_manual(values = c("#ffd700", "#d73027"))+
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~clusters, ncol=4)

plot(g3_3)


g_3_2_2 <- ggplot(e_3, aes(x=phase, y=mean_avpar, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=mean_avpar-sd_avpar, ymax=mean_avpar+sd_avpar), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Paramutable loci insertions")+
  xlab("Phase")+
  theme(legend.position="none")+
  scale_y_continuous(limits = c(-0.5, 12.5))+
  scale_fill_manual(values = c("#ffd700", "#d73027"))+
  facet_wrap(~sampleid, ncol=8,labeller = labeller(sampleid = 
                                                     c("p1" = "clu=1%",
                                                       "p3" = "clu=3%",
                                                       "p10" = "clu=10%",
                                                       "p50" = "clu=50%",
                                                       "p1_10" = "clu=1% para=10%",
                                                       "p3_10" = "clu=3% para=10%",
                                                       "p10_10" = "clu=10% para=10%",
                                                       "p50_10" = "clu=50% para=10%")))

g2_2 / g_3_2_2
