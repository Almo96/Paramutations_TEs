library(ggplot2)
library(plyr)

setwd("/Users/ascarpa/Paramutations_TEs/Validation/Raw")

validation<-read.table("2022_08_01_Validation_5_trigger_loci", fill = TRUE, sep = "\t")
names(validation)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
                     "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

data_new <- validation
data_new$sampleid <- factor(data_new$sampleid,
                            levels = c("pt1", "pt2", "pt3", "pt4", "pt5", "pt6", "pt7", "pt8", "pt9"))

gl<-ggplot()+geom_line(data=data_new,aes(x=gen,group=rep,y=avtes*1000),alpha=0.4)+scale_y_log10()+theme(legend.position="none")+ylab("TE copies in population")+xlab("generation")+facet_wrap(~sampleid, ncol=3)
plot(gl)