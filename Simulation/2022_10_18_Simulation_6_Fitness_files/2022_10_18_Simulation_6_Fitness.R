library(tidyverse)
library(ggplot2)
theme_set(theme_bw())

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")


df_0.02<-read.table("2022_10_18_Simulation_6_Fitness_0.02", fill = TRUE, sep = "\t")
df_0.03<-read.table("2022_10_18_Simulation_6_Fitness_0.03", fill = TRUE, sep = "\t")
df_0.05<-read.table("2022_10_18_Simulation_6_Fitness_0.05", fill = TRUE, sep = "\t")
df_0.02_noxclu<-read.table("2022_10_18_Simulation_6_Fitness_0.02_noxclu", fill = TRUE, sep = "\t")
df_0.03_noxclu<-read.table("2022_10_18_Simulation_6_Fitness_0.03_noxclu", fill = TRUE, sep = "\t")
df_0.05_noxclu<-read.table("2022_10_18_Simulation_6_Fitness_0.05_noxclu", fill = TRUE, sep = "\t")

naming <- c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed",
            "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
            "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid", "extra")
names(df_0.02) <- naming
names(df_0.03) <- naming
names(df_0.05) <- naming
names(df_0.02_noxclu) <- naming
names(df_0.03_noxclu) <- naming
names(df_0.05_noxclu) <- naming


g_A_0.02 <- ggplot(df_0.02,aes(x = gen, y = avtes))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.02" = "Paramutable loci = 0%",
                                                        "p1_x0.02" = "Paramutable loci = 1%",
                                                        "p10_x0.02" = "Paramutable loci = 10%",
                                                        "p100_x0.02" = "Paramutable loci = 100%")))
plot(g_A_0.02)

g_B_0.02 <- ggplot(df_0.02,aes(x = gen, y = avw))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.02" = "Paramutable loci = 0%",
                                                        "p1_x0.02" = "Paramutable loci = 1%",
                                                        "p10_x0.02" = "Paramutable loci = 10%",
                                                        "p100_x0.02" = "Paramutable loci = 100%")))
plot(g_B_0.02)


g_A_0.03 <- ggplot(df_0.03,aes(x = gen, y = avtes))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.03" = "Paramutable loci = 0%",
                                                        "p1_x0.03" = "Paramutable loci = 1%",
                                                        "p10_x0.03" = "Paramutable loci = 10%",
                                                        "p100_x0.03" = "Paramutable loci = 100%")))
plot(g_A_0.03)

g_B_0.03 <- ggplot(df_0.03,aes(x = gen, y = avw))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.03" = "Paramutable loci = 0%",
                                                        "p1_x0.03" = "Paramutable loci = 1%",
                                                        "p10_x0.03" = "Paramutable loci = 10%",
                                                        "p100_x0.03" = "Paramutable loci = 100%")))
plot(g_B_0.03)


g_A_0.05 <- ggplot(df_0.05,aes(x = gen, y = avtes))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.05" = "Paramutable loci = 0%",
                                                        "p1_x0.05" = "Paramutable loci = 1%",
                                                        "p10_x0.05" = "Paramutable loci = 10%",
                                                        "p100_x0.05" = "Paramutable loci = 100%")))
plot(g_A_0.05)

g_B_0.05 <- ggplot(df_0.05,aes(x = gen, y = avw))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.05" = "Paramutable loci = 0%",
                                                        "p1_x0.05" = "Paramutable loci = 1%",
                                                        "p10_x0.05" = "Paramutable loci = 10%",
                                                        "p100_x0.05" = "Paramutable loci = 100%")))
plot(g_B_0.05)



g_A_0.02_noxclu <- ggplot(df_0.02_noxclu,aes(x = gen, y = avtes))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.02" = "Paramutable loci = 0%",
                                                        "p1_x0.02" = "Paramutable loci = 1%",
                                                        "p10_x0.02" = "Paramutable loci = 10%",
                                                        "p100_x0.02" = "Paramutable loci = 100%")))
plot(g_A_0.02_noxclu)

g_B_0.02_noxclu <- ggplot(df_0.02_noxclu,aes(x = gen, y = avw))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.02" = "Paramutable loci = 0%",
                                                        "p1_x0.02" = "Paramutable loci = 1%",
                                                        "p10_x0.02" = "Paramutable loci = 10%",
                                                        "p100_x0.02" = "Paramutable loci = 100%")))
plot(g_B_0.02_noxclu)


g_A_0.03_noxclu <- ggplot(df_0.03_noxclu,aes(x = gen, y = avtes))+
  geom_line(alpha=0.1, size=0.05)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid = 
                                                      c("p0_x0.03" = "Paramutable loci = 0%",
                                                        "p1_x0.03" = "Paramutable loci = 1%",
                                                        "p10_x0.03" = "Paramutable loci = 10%",
                                                        "p100_x0.03" = "Paramutable loci = 100%")))
plot(g_A_0.03_noxclu)

g_B_0.03_noxclu <- ggplot(df_0.03_noxclu,aes(x = gen, y = avw))+
  geom_line(alpha=0.1, size=0.05)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.03" = "Paramutable loci = 0%",
                                                        "p1_x0.03" = "Paramutable loci = 1%",
                                                        "p10_x0.03" = "Paramutable loci = 10%",
                                                        "p100_x0.03" = "Paramutable loci = 100%")))
plot(g_B_0.03_noxclu)

df_summary <- df_0.03_noxclu %>% 
  group_by(sampleid, rep) %>%
  summarize(min_fitness = min(avw))

df_summary$sampleid[df_summary$sampleid == "p0_x0.03"] <- "0% (Trap model)"
df_summary$sampleid[df_summary$sampleid == "p1_x0.03"] <- "1%"
df_summary$sampleid[df_summary$sampleid == "p10_x0.03"] <- "10%"
df_summary$sampleid[df_summary$sampleid == "p100_x0.03"] <- "100%"

boxplot(df_summary$min_fitness ~ df_summary$sampleid, xlab="Paramutable loci", ylab = "Minimum fitness")


g_A_0.05_noxclu <- ggplot(df_0.05_noxclu,aes(x = gen, y = avtes))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.05" = "Paramutable loci = 0%",
                                                        "p1_x0.05" = "Paramutable loci = 1%",
                                                        "p10_x0.05" = "Paramutable loci = 10%",
                                                        "p100_x0.05" = "Paramutable loci = 100%")))
plot(g_A_0.05_noxclu)

g_B_0.05_noxclu <- ggplot(df_0.05_noxclu,aes(x = gen, y = avw))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.05" = "Paramutable loci = 0%",
                                                        "p1_x0.05" = "Paramutable loci = 1%",
                                                        "p10_x0.05" = "Paramutable loci = 10%",
                                                        "p100_x0.05" = "Paramutable loci = 100%")))
plot(g_B_0.05_noxclu)