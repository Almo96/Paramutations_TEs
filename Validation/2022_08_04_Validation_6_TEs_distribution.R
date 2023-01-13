library(ggplot2)
library(plyr)
library(patchwork)

setwd("/Users/ascarpa/Paramutations_TEs/Validation/Raw")


t_1<-read.table("validation_6_1_mhp", fill = TRUE, sep = "\t")
names(t_1)<-c("rep","gen","chr","pos","locus","popfreq")
t_1$rep<-as.factor(t_1$rep)
t_1$gen<-as.factor(t_1$gen)
t_1<-subset(t_1,rep==1)
t_1<-subset(t_1,gen==0 | gen==25 | gen==50 | gen==75 | gen==100)
g_1<-ggplot(data=t_1,aes(x=pos))+geom_histogram(binwidth=10000)+facet_grid(gen~chr, scales="free_x", space = "free_x")+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c("0","0.5m","1m","1.5m"))+
  xlab("position")+ylab("counter per 10kb bin log10")+scale_y_log10()
plot(g_1)


t_2<-read.table("validation_6_2_mhp", fill = TRUE, sep = "\t")
names(t_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_2$rep<-as.factor(t_2$rep)
t_2$gen<-as.factor(t_2$gen)
t_2<-subset(t_2,rep==1)
t_2<-subset(t_2,gen==0 | gen==25 | gen==50 | gen==75 | gen==100)
g_2<-ggplot(data=t_2,aes(x=pos, fill=locus))+geom_histogram(binwidth=10000)+facet_grid(gen~chr, scales="free_x", space = "free_x")+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c("0","0.5m","1m","1.5m"))+
  xlab("position")+ylab("counter per 10kb bin")
plot(g_2)

t_2_2<-read.table("validation_6_2_mhp", fill = TRUE, sep = "\t")
names(t_2_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_2_2$gen<-as.factor(t_2_2$gen)
t_2_2<-subset(t_2_2, gen==0)
g_2_2<-ggplot()+
  geom_bar(data=t_2_2,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(legend.position="none")
  ylab("relative frequencies")

t_2_3<-read.table("validation_6_2_mhp", fill = TRUE, sep = "\t")
names(t_2_3)<-c("rep","gen","chr","pos","locus","popfreq")
t_2_3$gen<-as.factor(t_2_3$gen)
t_2_3<-subset(t_2_3, gen==100)
g_2_3<-ggplot()+
  geom_bar(data=t_2_3,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("relative frequencies")

g_2_2+g_2_3


t_3<-read.table("validation_6_3_mhp", fill = TRUE, sep = "\t")
names(t_3)<-c("rep","gen","chr","pos","locus","popfreq")
t_3$rep<-as.factor(t_3$rep)
t_3$gen<-as.factor(t_3$gen)
t_3<-subset(t_3,rep==1)
t_3<-subset(t_3, gen==0 | gen==25 | gen==50 | gen==75 | gen==100)
g_3<-ggplot(data=t_3,aes(x=pos, fill=locus))+geom_histogram(binwidth=10000)+facet_grid(gen~chr, scales="free_x", space = "free_x")+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c("0","0.5m","1m","1.5m"))+
  xlab("position")+ylab("counter per 10kb bin")
plot(g_3)

t_3_2<-read.table("validation_6_3_mhp", fill = TRUE, sep = "\t")
names(t_3_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_3_2$gen<-as.factor(t_3_2$gen)
t_3_2<-subset(t_3_2, gen==0)
g_3_2<-ggplot()+
  geom_bar(data=t_3_2,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(legend.position="none")+
  ylab("relative frequencies")

t_3_3<-read.table("validation_6_3_mhp", fill = TRUE, sep = "\t")
names(t_3_3)<-c("rep","gen","chr","pos","locus","popfreq")
t_3_3$gen<-as.factor(t_3_3$gen)
t_3_3<-subset(t_3_3, gen==100)
g_3_3<-ggplot()+
  geom_bar(data=t_3_3,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("relative frequencies")

g_3_2+g_3_3

t_4<-read.table("validation_6_4_mhp", fill = TRUE, sep = "\t")
names(t_4)<-c("rep","gen","chr","pos","locus","popfreq")
t_4$rep<-as.factor(t_4$rep)
t_4$gen<-as.factor(t_4$gen)
t_4<-subset(t_4,rep==1)
t_4<-subset(t_4, gen==0 | gen==25 | gen==50 | gen==75 | gen==100)
g_4<-ggplot(data=t_4,aes(x=pos, fill=locus))+geom_histogram(binwidth=10000)+facet_grid(gen~chr, scales="free_x", space = "free_x")+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c("0","0.5m","1m","1.5m"))+
  xlab("position")+ylab("counter per 10kb bin")
plot(g_4)

t_4_2<-read.table("validation_6_4_mhp", fill = TRUE, sep = "\t")
names(t_4_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_4_2$gen<-as.factor(t_4_2$gen)
t_4_2<-subset(t_4_2, gen==0)
g_4_2<-ggplot()+
  geom_bar(data=t_4_2,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(legend.position="none")+
  ylab("relative frequencies")

t_4_3<-read.table("validation_6_4_mhp", fill = TRUE, sep = "\t")
names(t_4_3)<-c("rep","gen","chr","pos","locus","popfreq")
t_4_3$gen<-as.factor(t_4_3$gen)
t_4_3<-subset(t_4_3, gen==100)
g_4_3<-ggplot()+
  geom_bar(data=t_4_3,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("relative frequencies")

g_4_2+g_4_3


t_5<-read.table("validation_6_5_mhp", fill = TRUE, sep = "\t")
names(t_5)<-c("rep","gen","chr","pos","locus","popfreq")

t_5$rep<-as.factor(t_5$rep)
t_5$gen<-as.factor(t_5$gen)
t_5<-subset(t_5,rep==4)
t_5<-subset(t_5, gen==0 | gen==25 | gen==50 | gen==75 | gen==100)

g_5<-ggplot(data=t_5,aes(x=pos, fill=locus))+geom_histogram(binwidth=10000)+facet_grid(gen~chr, scales="free_x", space = "free_x")+
  scale_x_continuous(breaks=c(0,500000,1000000,1500000),labels=c("0","0.5m","1m","1.5m"))+
  xlab("position")+ylab("counter per 10kb bin")

plot(g_5)


t_5_2<-read.table("validation_6_5_mhp", fill = TRUE, sep = "\t")
names(t_5_2)<-c("rep","gen","chr","pos","locus","popfreq")
t_5_2$gen<-as.factor(t_5_2$gen)
t_5_2<-subset(t_5_2, gen==0)
g_5_2<-ggplot()+
  geom_bar(data=t_5_2,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  ylab("relative frequencies")

t_5_3<-read.table("validation_6_5_mhp", fill = TRUE, sep = "\t")
names(t_5_3)<-c("rep","gen","chr","pos","locus","popfreq")
t_5_3$gen<-as.factor(t_5_3$gen)
t_5_3<-subset(t_5_3, gen==100)
g_5_3<-ggplot()+
  geom_bar(data=t_5_3,aes(x = locus, y = (..count..)/sum(..count..), fill = locus))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(legend.position="none")+
  ylab("relative frequencies")

g_5_2+g_5_3
