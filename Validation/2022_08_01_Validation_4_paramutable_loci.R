library(ggplot2)
library(dplyr)
theme_set(theme_bw())

setwd("/Users/ascarpa/Paramutations_TEs/Validation/Raw")

validation<-read.table("2022_08_01_Validation_4_paramutable_loci", fill = TRUE, sep = "\t")
names(validation)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
                       "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")



data_new <- validation
data_new$sampleid <- factor(data_new$sampleid,
                            levels = c("ps0", "ps50", "ps100", "ps100_0", "ps100_50", "ps100_100"))

gl<-ggplot()+geom_line(data=data_new,aes(x=gen,group=rep,y=avtes*1000),alpha=0.4)+scale_y_log10()+theme(legend.position="none")+ylab("TE copies in population")+xlab("generation")+facet_wrap(~sampleid)
plot(gl)
