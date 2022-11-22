2022_11_18_Simulation_8\_Paramutation_fix
================
Almo
2022-11-18

## Introduction

With this simulation we wanted to study the possibility of fixation of a
paramutation in a population.

### Initial conditions:

A population of 1000, 1 chromosomes of size 10 Mb, no piRNA clusters and
a TE for each chromosome in position 1, the same position is also a
paramutable locus. Half of the population with maternal piRNA
deposition.

We used 1000 replicates for each simulation.

## Materials & Methods

version: invadego0.23

-   seed 8_1: 1668790104187535070

-   seed 8_2: 1669023324203668900

-   seed 8_3: 1669038874472705808

-   seed 8_4: 1669051103177746034

-   seed 8_5: 1669063981462188015

-   seed 8_6: 1669049996301417348

### Commands for the simulation:

``` bash
echo "500 R 1;1;1
500 R 0;1;1" > 2022_11_18_input_08

folder="/Users/ascarpa/Paramutations_TEs/Simulation/Raw"
tool="/Users/ascarpa/invade-invadego/invadego023"

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0 --gen 5000 --genome mb:10 --steps 5000 --rr 4 --paramutation 999999:1 --rep 1000 --silent > $folder/2022_11_18_simulation_8_1

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.01 -x 0.01 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_2

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.01 -x 0.1 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_3

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.1 -x 0.01 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_4

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.1 -x 0.1 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_5

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.1 -x 0.001 --gen 5000 --genome mb:10 --steps 5000 --rr 4 --paramutation 999999:1 --rep 1000 --silent > $folder/2022_11_18_simulation_8_6
```

### Visualization in R

Setting the environment

``` r
library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
theme_set(theme_bw())
```

Visualization:

``` r
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
  theme(legend.position="none")+
  theme_bw()

plot(g)
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
df_2 <- data.frame()
df_2 <- subset(df, gen == 5000)

df_2_1 <- df_2
df_2_1$fwpar_yespi <- c("No piRNAs", "Fixed paramutation")


g_2<-ggplot(df_2_1, aes(x=as.factor(fwpar_yespi), fill=as.factor(fwpar_yespi) )) + 
  geom_bar( ) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position="none")

plot(g_2)
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
df_3 <- df_2 %>%
  dplyr::count(fwpar_yespi)


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
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
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
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

``` r
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
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-5.png)<!-- -->

``` r
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
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-6.png)<!-- -->

``` r
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
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-7.png)<!-- -->

``` r
#--cluster kb:0 --u 0.1 -x 0.001
df_x1_5 <- read.table("2022_11_18_Simulation_8_6", fill = TRUE, sep = "\t")
names(df_x1_5) <- names_vector
df_x2_5 <- data.frame()
df_x2_5 <- subset(df_x1_3, gen == 5000)

df_x3_5 <- df_x2_5 %>%
  dplyr::count(fwpar_yespi)

df_x3_5$fwpar_yespi <- c("Fixed paramutation")


g_x5 <-ggplot(df_x3_5, aes(x="", y=n, fill=factor(fwpar_yespi)))+
  geom_col(color = "black") +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())+
  ggtitle("--cluster kb:0 --u 0.1 -x 0.001")

plot(g_x5)
```

![](2022_11_18_Simulation_8_Paramutation_fix_files/figure-gfm/unnamed-chunk-3-8.png)<!-- -->

## Conclusions

A paramutation can be fixed in the population under selection for TEs.
