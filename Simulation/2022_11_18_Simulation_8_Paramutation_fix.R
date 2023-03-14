library(tidyverse)
library(ggplot2)
library(dplyr)
library(plyr)
library(patchwork)
library(ggpubr)
theme_set(theme_bw())


setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

df<-read.table("2022_11_18_Simulation_8_1", fill = TRUE, sep = "\t")

names_vector <- c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
             "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
             "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6")

names(df) <- names_vector

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))

g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=fwpar_yespi,group=rep),alpha=1,size=0.7)+
  xlab("generation")+
  ylab("")+
  theme(legend.position="none")

plot(g)


df_2 <- data.frame()
df_2 <- subset(df, gen == 5000)

df_2_1 <- df_2
df_2_1$fwpar_yespi <- c("No piRNAs", "Fixed paramutation")

df_3 <- df_2 %>%
  dplyr::count(fwpar_yespi)

g_2<-ggplot(df_3, aes(x=fwpar_yespi, y=n/10, fill=factor(fwpar_yespi)))+
  geom_col(color = "black")+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ggtitle("x = 0, u = 0")+
  ylab("populations [%]")+
  xlab("")+
  theme(legend.position = "none",
      plot.title = element_text(size=24),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20),
      axis.title.x = element_text(size=24),
      axis.title.y = element_text(size=24),
      strip.text = element_text(size = 24))

plot(g_2)

df_3$fwpar_yespi <- c("no piRNAs", "Fixed paramutation")

g_3 <- ggplot(df_3, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

plot(g_3)


#--cluster kb:0 --u 0.01 -x 0.01
df_x1 <- read.table("2022_11_18_Simulation_8_2", fill = TRUE, sep = "\t")
names(df_x1) <- names_vector
df_x2 <- data.frame()
df_x2 <- subset(df_x1, gen == 5000)

df_x3 <- df_x2 %>%
  dplyr::count(fwpar_yespi)

df_x3$fwpar_yespi <- c("No piRNAs", "Some piRNAs", "Fixed paramutation")


g_x1 <-ggplot(df_x3, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
        ggtitle("--cluster kb:0 --u 0.01 -x 0.01")

plot(g_x1)


#--cluster kb:0 --u 0.01 -x 0.1
df_x1_2 <- read.table("2022_11_18_Simulation_8_3", fill = TRUE, sep = "\t")
names(df_x1_2) <- names_vector
df_x2_2 <- data.frame()
df_x2_2 <- subset(df_x1_2, gen == 5000)

df_x3_2 <- df_x2_2 %>%
  dplyr::count(fwpar_yespi)

df_x3_2$fwpar_yespi <- c("No piRNAs", "Fixed paramutation")


g_x2 <-ggplot(df_x3_2, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
        ggtitle("--cluster kb:0 --u 0.01 -x 0.1")

plot(g_x2)


#--cluster kb:0 --u 0.1 -x 0.01
df_x1_3 <- read.table("2022_11_18_Simulation_8_4", fill = TRUE, sep = "\t")
names(df_x1_3) <- names_vector
df_x2_3 <- data.frame()
df_x2_3 <- subset(df_x1_3, gen == 5000)

df_x3_3 <- df_x2_3 %>%
  dplyr::count(fwpar_yespi)

df_x3_3$fwpar_yespi <- c("Fixed paramutation")


g_x3 <-ggplot(df_x3_3, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
        ggtitle("--cluster kb:0 --u 0.1 -x 0.01")

plot(g_x3)


#--cluster kb:0 --u 0.1 -x 0.1
df_x1_4 <- read.table("2022_11_18_Simulation_8_5", fill = TRUE, sep = "\t")
names(df_x1_4) <- names_vector
df_x2_4 <- data.frame()
df_x2_4 <- subset(df_x1_4, gen == 5000)

df_x3_4 <- df_x2_4 %>%
  dplyr::count(fwpar_yespi)

df_x3_4$fwpar_yespi <- c("Fixed paramutation")


g_x4 <-ggplot(df_x3_4, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
        ggtitle("--cluster kb:0 --u 0.1 -x 0.1")

plot(g_x4)


#--cluster kb:0 --u 0 -x 0
df_x1_5 <- read.table("2022_11_18_Simulation_8_6", fill = TRUE, sep = "\t")
names(df_x1_5) <- names_vector
df_x2_5 <- data.frame()
df_x2_5 <- subset(df_x1_5, gen == 5000)

df_x3_5 <- df_x2_5 %>%
  dplyr::count(fwpar_yespi)

df_x3_5$fwpar_yespi <- c("lost", "fixed")


g_x5 <-ggplot(df_x3_5, aes(x=fwpar_yespi, y=n/10, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ggtitle("x = 0, u = 0")+
  ylab("populations [%]")+
  xlab("")+
  theme(legend.position = "none",
        plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))

plot(g_x5)


#--cluster kb:0,0,0,0,0 --u 0.01 -x 0.01
df_x1_6 <- read.table("2022_11_18_simulation_8_7", fill = TRUE, sep = "\t")
names(df_x1_6) <- names_vector
df_x2_6 <- data.frame()
df_x2_6 <- subset(df_x1_6, gen == 5000)

df_x3_6 <- df_x2_6 %>%
  dplyr::count(fwpar_yespi)

df_x3_6$fwpar_yespi <- c("lost", "fixed")


g_x6 <-ggplot(df_x3_6, aes(x=fwpar_yespi, y=n/10, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ggtitle("x = 0.01, u = 0.01")+
  ylab("populations [%]")+
  xlab("")+
  theme(legend.position = "none",
        plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))

plot(g_x6)

g_x6_2 <-ggplot(df_x3_6, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
  ggtitle("5 chromosomes, 1 paramutable locus per haploid genome, no clusters --u 0.01 -x 0.01")

plot(g_x6_2)

#--cluster kb:0,0,0,0,0 --u 0.01 -x 0.1
df_x1_7 <- read.table("2022_11_18_simulation_8_8", fill = TRUE, sep = "\t")
names(df_x1_7) <- names_vector
df_x2_7 <- data.frame()
df_x2_7 <- subset(df_x1_7, gen == 5000)

df_x3_7 <- df_x2_7 %>%
  dplyr::count(fwpar_yespi)

df_x3_7$fwpar_yespi <- c("lost", "fixed")

g_x7 <-ggplot(df_x3_7, aes(x=fwpar_yespi, y=n/10, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ggtitle("x = 0.1, u = 0.01")+
  ylab("populations [%]")+
  xlab("")+
  theme(legend.position = "none",
        plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))

plot(g_x7)

g_x7_2 <-ggplot(df_x3_7, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
  ggtitle("5 chromosomes, 1 paramutable locus per haploid genome, no clusters --u 0.01 -x 0.1")

plot(g_x7_2)

#--cluster kb:0,0,0,0,0 --u 0.01 -x 0.2
df_x1_8 <- read.table("2022_11_18_simulation_8_9", fill = TRUE, sep = "\t")
names(df_x1_8) <- names_vector
df_x2_8 <- data.frame()
df_x2_8 <- subset(df_x1_8, gen == 5000)

df_x3_8 <- df_x2_8 %>%
  dplyr::count(fwpar_yespi)

df_x3_8$fwpar_yespi <- c("lost", "fixed")

g_x8 <-ggplot(df_x3_8, aes(x=fwpar_yespi, y=n/10, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  ggtitle("x = 0.2, u = 0.01")+
  ylab("populations [%]")+
  xlab("")+
  theme(legend.position = "none",
        plot.title = element_text(size=24),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=24),
        strip.text = element_text(size = 24))

plot(g_x8)


g_x8_2 <-ggplot(df_x3_8, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
  ggtitle("5 chromosomes, 1 paramutable locus per haploid genome, no clusters --u 0.01 -x 0.2")

plot(g_x8_2)


g_final<-ggarrange(g_x5, g_x6, g_x7, g_x8,
          ncol = 4, nrow = 1, align = ("v"),
          labels = c("B", "C", "D", "E"), heights = c(2,2), widths = c(2,2)
)

g_final/g_final

pdf(file = "/Users/ascarpa/Paramutations_TEs/Pictures_paper/Figure_5/Figure_5.pdf", width = 10, height = 7.5)
ggarrange(g_x5, g_x6, g_x7, g_x8,
          ncol = 4, nrow = 1, align = ("v"),
          labels = c("B", "C", "D", "E"), heights = c(2,2), widths = c(2,2)
)
dev.off()


#g_x5 x = 0
yespi <- sum(df_x2_5$fwpar_yespi == 1)
nopi <- sum(df_x2_5$fwpar_nopi == 1)
m = matrix(c(yespi, nopi,500, 500),nrow=2,byrow=T)
chisq.test(m)
#Pearson's Chi-squared test with Yates' continuity correction
#data:  m
#X-squared = 0.072004, df = 1, p-value = 0.7884

#g_x6 x = 0.01
yespi <- sum(df_x2_6$fwpar_yespi == 1)
nopi <- sum(df_x2_6$fwpar_nopi == 1)
m = matrix(c(yespi, nopi,500, 500),nrow=2,byrow=T)
chisq.test(m)
#Pearson's Chi-squared test with Yates' continuity correction
#data:  m
#X-squared = 5.6344, df = 1, p-value = 0.01761

#g_x6 x = 0.1
yespi <- sum(df_x2_7$fwpar_yespi == 1)
nopi <- sum(df_x2_7$fwpar_nopi == 1)
m = matrix(c(yespi, nopi,500, 500),nrow=2,byrow=T)
chisq.test(m)
#Pearson's Chi-squared test with Yates' continuity correction
#data:  m
#X-squared = 125.57, df = 1, p-value < 2.2e-16

#g_x6 x = 0.2
yespi <- sum(df_x2_8$fwpar_yespi == 1)
nopi <- sum(df_x2_8$fwpar_nopi == 1)
m = matrix(c(yespi, nopi,500, 500),nrow=2,byrow=T)
chisq.test(m)
#Pearson's Chi-squared test with Yates' continuity correction
#data:  m
#X-squared = 431.29, df = 1, p-value < 2.2e-16
