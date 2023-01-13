library(ggplot2)
library(RColorBrewer)
library(plyr)
library(patchwork)

setwd("/Users/ascarpa/Paramutations_TEs/Validation/Raw")


t_1<-read.table("validation_7_1_mhp", fill = TRUE, sep = "\t")
names(t_1)<-c("rep","gen","chr","pos","locus","popfreq")
t_1$rep<-as.factor(t_1$rep)
t_1$gen<-as.factor(t_1$gen)
t_1<-subset(t_1,rep==10)
t_1<-subset(t_1,gen==0 | gen==25 | gen==50)
g_1<-ggplot(data=t_1,aes(x=pos, fill=locus))+geom_histogram(binwidth=10000)+facet_grid(gen~chr, scales="free_x", space = "free_x")+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c("0","0.5m","1m","1.5m"))+
  xlab("position")+ylab("counter per 10kb bin")
plot(g_1)

t_1_2<-read.table("validation_7_1_mhp", fill = TRUE, sep = "\t")
names(t_1_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_1_2$gen<-as.factor(t_1_2$gen)
t_1_2<-subset(t_1_2, gen==0)
g_1_2<-ggplot()+
  geom_bar(data=t_1_2,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ggtitle("Generation 0")+
  theme(legend.position="none", plot.title = element_text(size=14, face="bold.italic"))+
ylab("relative frequencies")

t_1_3<-read.table("validation_7_1_mhp", fill = TRUE, sep = "\t")
names(t_1_3)<-c("rep","gen","chr","pos","locus","popfreq")
t_1_3$gen<-as.factor(t_1_3$gen)
t_1_3<-subset(t_1_3, gen==50)
g_1_3<-ggplot()+
  geom_bar(data=t_1_3,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ggtitle("Generation 50")+
  theme(plot.title = element_text(size=14, face="bold.italic"))+
  ylab("relative frequencies")

g_1_2+g_1_3


t_2<-read.table("validation_7_2_mhp", fill = TRUE, sep = "\t")
names(t_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_2$rep<-as.factor(t_2$rep)
t_2$gen<-as.factor(t_2$gen)
t_2<-subset(t_2,rep==1)
t_2<-subset(t_2,gen==0 | gen==25 | gen==50 | gen==75 | gen==100)
g_2<-ggplot(data=t_2,aes(x=pos, fill=locus))+geom_histogram(binwidth=10000)+facet_grid(gen~chr, scales="free_x", space = "free_x")+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c("0","0.5m","1m","1.5m"))+
  xlab("position")+ylab("counter per 10kb bin")
plot(g_2)

t_2_2<-read.table("validation_7_2_mhp", fill = TRUE, sep = "\t")
names(t_2_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_2_2$gen<-as.factor(t_2_2$gen)
t_2_2<-subset(t_2_2, gen==0)
g_2_2<-ggplot()+
  geom_bar(data=t_2_2,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ggtitle("Generation 0")+
  theme(legend.position="none", plot.title = element_text(size=14, face="bold.italic"))+
  ylab("relative frequencies")

t_2_3<-read.table("validation_7_2_mhp", fill = TRUE, sep = "\t")
names(t_2_3)<-c("rep","gen","chr","pos","locus","popfreq")
t_2_3$gen<-as.factor(t_2_3$gen)
t_2_3<-subset(t_2_3, gen==100)
g_2_3<-ggplot()+
  geom_bar(data=t_2_3,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ggtitle("Generation 100")+
  theme(plot.title = element_text(size=14, face="bold.italic"))+
  ylab("relative frequencies")

g_2_2+g_2_3


df_1<-read.table("2022_08_05_Validation_7_selection_1", fill = TRUE, sep = "\t")
names(df_1)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
              "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

gt_1<-ggplot()+
  theme(legend.position="none")+
  geom_line(data=df_1,aes(x=gen, y=avtes , group=rep,color=phase),alpha=1,size=0.7)+
  ylab("TE population frequency") + xlab("generation")+
  scale_y_log10()+
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                       c("psel3" = "x = 0",
                                         "psel4" = "x = 0.1",
                                         "psel5" = "x = 0.01",
                                         "psel6" = "x = 0.001",
                                         "psel7" = "x = 0.0001")))
plot(gt_1)

df1_s<-subset(df_1, sampleid=="psel3")
df2_s<-subset(df_1, sampleid=="psel4")
df3_s<-subset(df_1, sampleid=="psel5")
df4_s<-subset(df_1, sampleid=="psel6")
df5_s<-subset(df_1, sampleid=="psel7")

p0=0.5
p<-p0
q<-1-p0
t=1000
traj<- matrix(NA, ncol=3, nrow=max(t)+1)
traj[1,1]<-p
traj[1,2]<-0
s=0
g=1
wAA=1-(2*s)
wAa=1-s
waa=1
while(g<=max(t)){
  w<- (p^2)*wAA + (2*p*q*wAa) + (q^2)*waa
  p<-( p*(p*wAA + q*wAa) ) / w
  q<-1-p
  traj[g+1,1]<-p
  traj[g+1,2]<-g
  g<-g+1
}
traj<- as.data.frame(traj)
traj[,3]<- "black"
colnames(traj)<- c("freq", "generations", "color")

g_s_1<-ggplot()+ 
  geom_line(df1_s, mapping=aes(x=gen, y=avtes, group=rep), color="red")+
  geom_line(traj, mapping=aes(x=generations, y=(freq*2)), color="black")+
  labs(x="generation", y="frequency of TEs in the population")
plot(g_s_1)

p0=0.5
p<-p0
q<-1-p0
t=1000
traj<- matrix(NA, ncol=3, nrow=max(t)+1)
traj[1,1]<-p
traj[1,2]<-0
s=0.1
g=1
wAA=1-(2*s)
wAa=1-s
waa=1
while(g<=max(t)){
  w<- (p^2)*wAA + (2*p*q*wAa) + (q^2)*waa
  p<-( p*(p*wAA + q*wAa) ) / w
  q<-1-p
  traj[g+1,1]<-p
  traj[g+1,2]<-g
  g<-g+1
}
traj<- as.data.frame(traj)
traj[,3]<- "black"
colnames(traj)<- c("freq", "generations", "color")
g_s_2<-ggplot()+ 
  geom_line(df2_s, mapping=aes(x=gen, y=avtes, group=rep), color="red")+
  geom_line(traj, mapping=aes(x=generations, y=(freq*2)), color="black")+
  labs(x="generation", y="frequency of TEs in the population")
plot(g_s_2)

p0=0.5
p<-p0
q<-1-p0
t=1000
traj<- matrix(NA, ncol=3, nrow=max(t)+1)
traj[1,1]<-p
traj[1,2]<-0
s=0.01
g=1
wAA=1-(2*s)
wAa=1-s
waa=1
while(g<=max(t)){
  w<- (p^2)*wAA + (2*p*q*wAa) + (q^2)*waa
  p<-( p*(p*wAA + q*wAa) ) / w
  q<-1-p
  traj[g+1,1]<-p
  traj[g+1,2]<-g
  g<-g+1
}
traj<- as.data.frame(traj)
traj[,3]<- "black"
colnames(traj)<- c("freq", "generations", "color")
g_s_3<-ggplot()+ 
  geom_line(df3_s, mapping=aes(x=gen, y=avtes, group=rep), color="red")+
  geom_line(traj, mapping=aes(x=generations, y=(freq*2)), color="black")+
  labs(x="generation", y="frequency of TEs in the population")
plot(g_s_3)

p0=0.5
p<-p0
q<-1-p0
t=1000
traj<- matrix(NA, ncol=3, nrow=max(t)+1)
traj[1,1]<-p
traj[1,2]<-0
s=0.001
g=1
wAA=1-(2*s)
wAa=1-s
waa=1
while(g<=max(t)){
  w<- (p^2)*wAA + (2*p*q*wAa) + (q^2)*waa
  p<-( p*(p*wAA + q*wAa) ) / w
  q<-1-p
  traj[g+1,1]<-p
  traj[g+1,2]<-g
  g<-g+1
}
traj<- as.data.frame(traj)
traj[,3]<- "black"
colnames(traj)<- c("freq", "generations", "color")
g_s_4<-ggplot()+ 
  geom_line(df4_s, mapping=aes(x=gen, y=avtes, group=rep), color="red")+
  geom_line(traj, mapping=aes(x=generations, y=(freq*2)), color="black")+
  labs(x="generation", y="frequency of TEs in the population")
plot(g_s_4)

p0=0.5
p<-p0
q<-1-p0
t=1000
traj<- matrix(NA, ncol=3, nrow=max(t)+1)
traj[1,1]<-p
traj[1,2]<-0
s=0.0001
g=1
wAA=1-(2*s)
wAa=1-s
waa=1
while(g<=max(t)){
  w<- (p^2)*wAA + (2*p*q*wAa) + (q^2)*waa
  p<-( p*(p*wAA + q*wAa) ) / w
  q<-1-p
  traj[g+1,1]<-p
  traj[g+1,2]<-g
  g<-g+1
}
traj<- as.data.frame(traj)
traj[,3]<- "black"
colnames(traj)<- c("freq", "generations", "color")
g_s_5<-ggplot()+ 
  geom_line(df5_s, mapping=aes(x=gen, y=avtes, group=rep), color="red")+
  geom_line(traj, mapping=aes(x=generations, y=(freq*2)), color="black")+
  labs(x="generation", y="frequency of TEs in the population")
plot(g_s_5)


df_2<-read.table("2022_08_05_Validation_7_selection_2", fill = TRUE, sep = "\t")
names(df_2)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
            "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")
df_2$sampleid <- factor(df_2$sampleid, levels=c("psel8", "psel9", "psel10", "psel11","psel12"))

gt_2<-ggplot()+
  theme(legend.position="none")+
  geom_line(data=df_2, aes(x=gen, y=avtes , group=rep,color=phase),alpha=1,size=0.7)+
  ylab("TE population frequency") + xlab("generation")+
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("psel8" = "N = 10",
                                                "psel9" = "N = 100",
                                                "psel10" = "N = 1000",
                                                "psel11" = "N = 10000",
                                                "psel12" = "N = 100000")))
plot(gt_2)

gt_2_2<-ggplot()+
  theme(legend.position="none")+
  geom_line(data=df_2, aes(x=gen, y=fixed , group=rep,color=phase),alpha=1,size=0.7)+
  ylab("Fixed TEs") + xlab("generation")+
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("psel8" = "N = 10",
                                                "psel9" = "N = 100",
                                                "psel10" = "N = 1000",
                                                "psel11" = "N = 10000",
                                                "psel12" = "N = 100000")))
plot(gt_2_2)
