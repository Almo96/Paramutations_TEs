library(tidyverse)
library(ggplot2)
library(patchwork)
library(plotrix)
library(RColorBrewer)
library(ggpubr)
library(dplyr)


p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

df<-read.table("2022_08_09_Simulation_1_Paramutations", fill = TRUE, sep = "\t")
names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed",
             "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
             "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p0", "p1", "p10","p100"))


g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  theme( axis.text = element_text( size = 24 ),
               axis.text.x = element_text( size = 24 ),
               axis.title = element_text( size = 24, face = "bold" ),
               legend.position="none",
               # The new stuff
               strip.text = element_text(size = 24))+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                       c("p0" = "Paramutable loci = 0% (Trap model)",
                                         "p1" = "Paramutable loci = 1%",
                                         "p10" = "Paramutable loci = 10%",
                                         "p100" = "Paramutable loci = 100%")))
                                       
plot(g)

pdf(file = "/Users/ascarpa/Paramutations_TEs/Pictures_paper/Figure_1/Figure_1.pdf", width = 10, height = 7.5)
g
dev.off()

clus_ins <- aggregate(x = df$fwcli,
            by = list(df$gen, df$sampleid),
            FUN = sum)
names(clus_ins) <- c("gen", "sampleid", "percentcli")
df2<-subset(clus_ins, gen == 0 | gen == 100 | gen == 1000 | gen == 2500 | gen == 5000)

g2 <- ggplot(df2, aes(x=as.character(gen), y=percentcli)) + 
      geom_bar(stat = "identity", aes(fill=gen)) +
      ylab("Percentage of individuals with a cluster insertion")+
      xlab("generation")+
      theme(legend.position = "none")+
      theme_bw()+
      facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("p0" = "Paramutable loci = 0% (Trap model)",
                                                "p1" = "Paramutable loci = 1%",
                                                "p10" = "Paramutable loci = 10%",
                                                "p100" = "Paramutable loci = 100%")))

plot(g2)


setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")
df3<-read.table("2022_08_19_100_simulations", fill = TRUE, sep = "\t")
names(df3)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed",
             "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
             "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid", "extra")


df3 = subset(df3, gen == 5000 )
df3 <- df3 %>% 
  select(-c("extra"))

df3_2 <- df3 %>% 
  dplyr::group_by(sampleid) %>% 
  dplyr::summarize(av_cli = mean(avcli), sd_cli = sd(avcli), cv_cli_percent = sd(avcli)/mean(avcli), 
            av_tes = mean(avtes), sd_tes = sd(avtes), cv_tes_percent = sd(avtes)/mean(avtes))

percent_para<-c(0,10,100,15,20,25,30,35,40,45,5,50,55,60,65,70,75,80,85,90,95)
df3_2$sampleid<-percent_para
df3_2 <- df3_2[order(df3_2$sampleid),]


coeff_3=32.5
g3 <- ggplot(df3_2, aes(x=sampleid/100))+
  geom_point(aes(y=av_cli*coeff_3), color="blue")+
  geom_line(aes(y=av_cli*coeff_3), color="blue")+
  geom_point(aes(y=av_tes), color="red")+
  geom_line(aes(y=av_tes), color="red")+
  scale_x_continuous(labels = scales::percent)+
  ggtitle("3% piRNA clusters")+
  scale_y_continuous(
    name = "TEs insertions per diploid individual",
    sec.axis = sec_axis(~./coeff_3, name="cluster insertions per diploid individual")
  )+
  xlab("Paramutable loci")+
  theme(legend.position="none",
        plot.title = element_text(size=14, face="bold"),
        axis.title.y = element_text(color = "red", size=10),
        axis.title.y.right = element_text(color = "blue", size=10)
  )+
  theme_bw()

plot(g3)

g_3_2 <- ggplot(df3_2, aes(x=sampleid/100, y=cv_tes_percent))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se = FALSE)+
  stat_regline_equation(label.y = 0.18, aes(label = ..eq.label..))+
  stat_regline_equation(label.y = 0.16, aes(label = ..rr.label..))+
  xlab("Paramutable loci")+
  ylab("Coefficient of variation TEs insertions per individual")+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.4))+
  ggtitle("Effects of paramutations on the cv of TE insertions per individual")+
  theme(legend.position="none",
        plot.title = element_text(size=14, face="bold"))+
  theme_bw()

plot(g_3_2)

g_3_3 <- ggplot(df3_2, aes(x=sampleid/100, y=cv_cli_percent))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se = FALSE)+
  stat_regline_equation(label.y = 2.8, aes(label = ..eq.label..))+
  stat_regline_equation(label.y = 2.6, aes(label = ..rr.label..))+
  xlab("Paramutable loci")+
  ylab("Coefficient of variation cluster insertions per individual")+
  scale_x_continuous(labels = scales::percent)+
  scale_y_continuous(labels = scales::percent)+
  ggtitle("Effects of paramutations on the cv of cluster insertions")+
  theme(legend.position="none",
        plot.title = element_text(size=14, face="bold"))+
  theme_bw()

plot(g_3_3)


df3_order <- df3
df3_order <- df3 %>%
  mutate(sampleid = str_replace(sampleid,"para0", "0")) %>%
  mutate(sampleid = str_replace(sampleid,"para5", "5")) %>%
  mutate(sampleid = str_replace(sampleid,"para10", "10")) %>%
  mutate(sampleid = str_replace(sampleid,"para15", "15")) %>%
  mutate(sampleid = str_replace(sampleid,"para20", "20")) %>%
  mutate(sampleid = str_replace(sampleid,"para25", "25")) %>%
  mutate(sampleid = str_replace(sampleid,"para30", "30")) %>%
  mutate(sampleid = str_replace(sampleid,"para35", "35")) %>%
  mutate(sampleid = str_replace(sampleid,"para40", "40")) %>%
  mutate(sampleid = str_replace(sampleid,"para45", "45")) %>%
  mutate(sampleid = str_replace(sampleid,"para50", "50")) %>%
  mutate(sampleid = str_replace(sampleid,"para55", "55")) %>%
  mutate(sampleid = str_replace(sampleid,"para60", "60")) %>%
  mutate(sampleid = str_replace(sampleid,"para65", "65")) %>%
  mutate(sampleid = str_replace(sampleid,"para70", "70")) %>%
  mutate(sampleid = str_replace(sampleid,"para75", "75")) %>%
  mutate(sampleid = str_replace(sampleid,"para80", "80")) %>%
  mutate(sampleid = str_replace(sampleid,"para85", "85")) %>%
  mutate(sampleid = str_replace(sampleid,"para90", "90")) %>%
  mutate(sampleid = str_replace(sampleid,"para95", "95")) %>%
  mutate(sampleid = str_replace(sampleid,"para100", "100"))

df3_order$sampleid<-as.integer(df3_order$sampleid)
df3_order <- df3_order[order(df3_order$sampleid),]

g_bar_av_TEs <- ggplot(df3_order, aes(x=as.factor(sampleid), y=avtes)) + 
  geom_boxplot() +
  xlab("Size of paramutable loci as [%]") +
  ylab("TEs insertions per diploid individual")

plot(g_bar_av_TEs)

g_bar_clusters <- ggplot(df3_order, aes(x=as.factor(sampleid), y=avcli)) + 
  geom_boxplot() +
  xlab("Size of paramutable loci as [%]") +
  ylab("cluster insertions per diploid individual")

plot(g_bar_clusters)
