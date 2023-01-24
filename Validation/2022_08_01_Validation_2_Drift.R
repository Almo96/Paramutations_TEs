library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/Users/ascarpa/Paramutations_TEs/Validation/Raw")

validation<-read.table("2022_08_01_Validation_2_Drift", fill = TRUE, sep = "\t")
names(validation)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi","avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

data_1 <- validation[which(validation$sampleid == "pd250"),names(validation) %in% c("rep","fixed")]
gl_1<-ggplot(data_1,aes(x=fixed))+
      geom_histogram(binwidth = 1)+
      ylim(0,65)+
      ggtitle("Population 250")+
      geom_vline(xintercept=20, lwd=1,lty=2, colour="blue")+
      theme(plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))

data_2 <- validation[which(validation$sampleid == "pd500"),names(validation) %in% c("rep","fixed")]
gl_2<-ggplot(data_2,aes(x=fixed))+
      geom_histogram(binwidth = 1)+
      ylim(0,65)+
      ggtitle("Population 500")+
      geom_vline(xintercept=10, lwd=1,lty=2, colour="blue")+
      theme(plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))

data_3 <- validation[which(validation$sampleid == "pd1000"),names(validation) %in% c("rep","fixed")]
gl_3<-ggplot(data_3,aes(x=fixed))+
      geom_histogram(binwidth = 1)+
      ylim(0,65)+
      ggtitle("Population 1000")+
      geom_vline(xintercept=5, lwd=1,lty=2, colour="blue")+
      theme(plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))


gl_1+gl_2+gl_3
