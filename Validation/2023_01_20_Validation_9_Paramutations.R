library(ggplot2)
library(plyr)
theme_set(theme_bw())


setwd("/Users/ascarpa/Paramutations_TEs/Validation/Raw")

df_1<-read.table("validation_9_1", fill = TRUE, sep = "\t")
df_2<-read.table("validation_9_2", fill = TRUE, sep = "\t")
df_3<-read.table("validation_9_3", fill = TRUE, sep = "\t")


naming<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
             "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
             "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6")

names(df_1)<-naming
names(df_2)<-naming
names(df_3)<-naming

df_1$sampleid <- "p0"
df_2$sampleid <- "p50"
df_3$sampleid <- "p100"

df_total <- rbind(df_1, df_2, df_3)

df_total$sampleid <- factor(df_total$sampleid,
                            levels = c("p0", "p50", "p100"))

g<-ggplot()+geom_line(data=df_total,aes(x=gen,group=rep,y=avtes*1000),alpha=0.4)+
  theme(legend.position="none")+
  ylab("TEs insertions per diploid individual")+
  xlab("generation")+
  theme(plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))+
  facet_wrap(~sampleid, ncol=3, labeller = labeller(sampleid = 
                                                      c("p0" = "p = 0%",
                                                        "p50" = "p = 50%",
                                                        "p100" = "p = 100%")))

plot(g)
