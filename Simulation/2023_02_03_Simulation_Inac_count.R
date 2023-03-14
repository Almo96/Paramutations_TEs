library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")

df <- read.table("2023_03_02_Simulation_supp_inac", fill = TRUE, sep = "\t")
names(df) <- c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
               "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
               "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p0", "p1", "p10", "p100", "p30_0.1", "p30_1", "p30_10", "p30_100"))

df <- df[,-29]

count_5000 <- df %>%
  filter(gen == 5000) %>%
  group_by(sampleid) %>%
  summarize(count_n = sum(phase == "inac"), total = n()) %>%
  ungroup() %>%
  mutate(ratio = count_n / total)

count_10000 <- df %>%
  filter(gen == 10000) %>%
  group_by(sampleid) %>%
  summarize(count_n = sum(phase == "inac"), total = n()) %>%
  ungroup() %>%
  mutate(ratio = count_n / total)

count_20000 <- df %>%
  filter(gen == 20000) %>%
  group_by(sampleid) %>%
  summarize(count_n = sum(phase == "inac"), total = n()) %>%
  ungroup() %>%
  mutate(ratio = count_n / total)
