library(tidyverse)
library(RColorBrewer)
library(ggpubr)
theme_set(theme_bw())

p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

df<-read.table("2023_05_12_Supplement_recombination", fill = TRUE, sep = "\t")
names(df) <- c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
               "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
               "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p1", "p2", "p3","p4"))


g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  theme( axis.text = element_text(size = 24),
         axis.text.x = element_text(size = 24),
         axis.title = element_text(size = 24),
         legend.position="none",
         strip.text = element_text(size = 24))+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("p1" = "r = 0",
                                                "p2" = "r = 4",
                                                "p3" = "r = 20",
                                                "p4" = "r = 40")))

plot(g)




df1 <- subset(df, phase %in% c("shot", "inac"))
df2 <- data.frame()

#new dataframe with only the first shotgun & the first inactive phase of each replicate
repcheck = 1
x = 1
y = 1
while (x<nrow(df1)) {
  if (repcheck != df1[x, 1]){
    y = 1
  }
  if (y == 1){
    if(df1[x, 13]  == "shot"){
      df2<-rbind(df2,df1[x,])
      y = 2
      repcheck = df1[x, 1]
    }
  }
  if (y == 2){
    if(df1[x, 13] == "inac"){
      df2<-rbind(df2,df1[x,])
      y = 1
    }
  }
  x = x+1
}

#Summary statistics
df2<-select (df2,-c(29))

df_count <- df2 %>%
  dplyr::count(sampleid, phase)

df_summary <- df2 %>% 
  dplyr::group_by(sampleid, phase) %>%
  dplyr::summarize(av_fwcli = mean(fwcli), sd_fwcli = sd(fwcli),
                   av_cli = mean(avcli), sd_cli = sd(avcli), cv_cli_percent = sd(avcli)/mean(avcli),
                   av_tes = mean(avtes), sd_tes = sd(avtes), cv_tes_percent = sd(avtes)/mean(avtes),
                   av_par = mean(avpar), sd_par = sd(avpar),
                   av_fwpar_yespi = mean(fwpar_yespi), sd_fwpar_yespi = sd(fwpar_yespi),
                   length_previous_phase = mean(gen),
                   sd_gen_phases = sd(gen))

df_summary <- cbind(df_count$n, df_summary)
colnames(df_summary)[1] ="n"

#CI 95%: z* sd/sqrt(population)
df_summary$ci_fwcli <- qt(0.975,df=df_summary$n-1)*(df_summary$sd_fwcli/sqrt(df_summary$n))
df_summary$ci_cli <- qt(0.975,df=df_summary$n-1)*(df_summary$sd_cli/sqrt(df_summary$n))
df_summary$ci_tes <- qt(0.975,df=df_summary$n-1)*(df_summary$sd_tes/sqrt(df_summary$n))
df_summary$ci_par <- qt(0.975,df=df_summary$n-1)*(df_summary$sd_par/sqrt(df_summary$n))
df_summary$ci_fwpar_yespi <- qt(0.975,df=df_summary$n-1)*(df_summary$sd_fwpar_yespi/sqrt(df_summary$n))

g_c <- ggplot(df_summary, aes(x=phase, y=av_tes, fill = phase)) + 
  geom_bar(stat = "identity") +
  geom_errorbar( aes(x=phase, ymin=av_tes-sd_tes, ymax=av_tes+sd_tes), width=0.2, colour="black", alpha=0.9, size=0.8)+
  ylab("Number of TE insertions per individual")+
  xlab("Phase")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)))+
  scale_fill_manual(values = c("#ffd700", "#d73027"))+
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid = 
                                                      c("p1" = "r = 0",
                                                        "p2" = "r = 4",
                                                        "p3" = "r = 20",
                                                        "p4" = "r = 40")))

plot(g_c)


df_5000 <- subset(df, gen ==5000)

g_5000 <- ggplot(df_5000, aes(x = sampleid, y = avtes)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("p1", "p2")), map_signif_level = TRUE, textsize = 2) +
  geom_signif(comparisons = list(c("p2", "p3")), map_signif_level = TRUE, textsize = 2) +
  geom_signif(comparisons = list(c("p3", "p4")), map_signif_level = TRUE, textsize = 2) +
  labs(x = "year", y = "copy number")

plot(g_5000)




#Removal of 5 clusters
df_new<-read.table("2022_12_10_Simulation_10_Clusters_removal_3000_new", fill = TRUE, sep = "\t")
names(df_new)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
                 "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
                 "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")


df_new$phase <- factor(df_new$phase, levels=c("rapi", "trig", "shot", "inac"))
df_new$sampleid <- factor(df_new$sampleid, levels=c("p0_0", "p0_1", "p0_2", "p0_3", "p0_4", "p0_5", "p10_0", "p10_1", "p10_2", "p10_3", "p10_4", "p10_5"))

g_new<-ggplot()+
  geom_line(data=df_new,aes(x=gen,y=avtes,group=rep,color=phase), alpha = 1, linewidth = 0.7)+
  geom_vline(xintercept = 3000, linetype="dashed", color = "black", linewidth = 0.7)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  scale_x_continuous(breaks = seq(0, 5000, by = 2500))+
  facet_wrap(~sampleid, ncol = 6, nrow = 2, labeller = labeller(sampleid = 
                                                                  c("p0_0" = "para = 0, clu -0 ",
                                                                    "p0_1" = "para = 0, clu -1 ",
                                                                    "p0_2" = "para = 0, clu -2 ",
                                                                    "p0_3" = "para = 0, clu -3 ",
                                                                    "p0_4" = "para = 0, clu -4 ",
                                                                    "p0_5" = "para = 0, clu -5 ",
                                                                    "p10_0" = "para = 10, clu -0 ",
                                                                    "p10_1" = "para = 10, clu -1 ",
                                                                    "p10_2" = "para = 10, clu -2 ",
                                                                    "p10_3" = "para = 10, clu -3 ",
                                                                    "p10_4" = "para = 10, clu -4 ",
                                                                    "p10_5" = "para = 10, clu -5 ")))


plot(g_new)

df_0_5 <- subset(df_new, sampleid == "p0_5" & gen == 5000)
df_10_5 <- subset(df_new, sampleid == "p10_5" & gen == 5000)

nrow(df_0_5)
nrow(df_10_5)
