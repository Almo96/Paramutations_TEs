2022_08_09_Simulation_1\_Paramutations
================
Almorò Scarpa

## Introduction

In this simulation we compared the effects of different shares of
trigger loci in a genome made of:

-   30% paramutable loci & 0% piRNA clusters,

-   10% paramutable loci & 0% piRNA clusters,

-   10% paramutable loci & 3% piRNA clusters.

### Initial conditions:

A population of 1000, 5 chromosomes of size 10 Mb, 30/10% of the genome
as paramutable loci and an initial number of TEs in the population equal
to 100.

We used 300 replicates for each simulation.

## Materials & Methods:

version: invadego0.2.1

-   seed p30_0.1: 1660131023489346000

-   seed p30_1: 1660120153070897000

-   seed p30_10: 1660120153071664000

-   seed p30_100: 1660120153072163000

-   seed p10_0.1: 1673959003642697000

-   seed p10_1: 1673984685519535000

-   seed p10_10: 1673995071515068000

-   seed p10_100: 1674003253780198000

-   seed p3_10_0.1: 1674011129204327000

-   seed p3_10_1: 1674018839215006000

-   seed p3_10_10: 1674026215672953000

-   seed p3_10_100: 1674032787174268000

### Commands for the simulation:

``` bash
folder="/Users/ascarpa/Paramutations_TEs/Simulation/Raw"
tool="/Users/ascarpa/invade-invadego/invadego022"

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1,3,5 --trigger 1000:1 --steps 20 --sampleid p30_0.1 > $folder/2022_08_09_simulation_2_1

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1,3,5 --trigger 100:1 --steps 20 --sampleid p30_1 > $folder/2022_08_09_simulation_2_2

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1,3,5 --trigger 10:1 --steps 20 --sampleid p30_10 > $folder/2022_08_09_simulation_2_3

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1,3,5 --trigger 1:0 --steps 20 --sampleid p30_100 > $folder/2022_08_09_simulation_2_4

cat 2022_08_09_simulation_2_1 2022_08_09_simulation_2_2 2022_08_09_simulation_2_3 2022_08_09_simulation_2_4 |grep -v "^Invade"|grep -v "^#" > 2022_08_09_Simulation_2_Trigger

folder="/Users/ascarpa/Paramutations_TEs/Simulation/Raw"
tool="/Users/ascarpa/invade-invadego/invadego023"

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 1000:1 --steps 20 --sampleid p10_0.1 > $folder/2023_01_17_simulation_2_5

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 100:1 --steps 20 --sampleid p10_1 > $folder/2023_01_17_simulation_2_6

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 10:1 --steps 20 --sampleid p10_10 > $folder/2023_01_17_simulation_2_7

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 1:0 --steps 20 --sampleid p10_100 > $folder/2023_01_17_simulation_2_8

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 1000:1 --steps 20 --sampleid p3_10_0.1 > $folder/2023_01_17_simulation_2_9

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 100:1 --steps 20 --sampleid p3_10_1 > $folder/2023_01_17_simulation_2_10

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 10:1 --steps 20 --sampleid p3_10_10 > $folder/2023_01_17_simulation_2_11

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --trigger 1:0 --steps 20 --sampleid p3_10_100 > $folder/2023_01_17_simulation_2_12
```

### Visualization in R

Setting the environment

``` r
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
theme_set(theme_bw())
```

Visualization: comparing the simulations with the prediction

``` r
p<-c("grey","#1a9850","#ffd700","#d73027")

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")
df<-read.table("2022_08_09_Simulation_2_Trigger", fill = TRUE, sep = "\t")

names(df)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
             "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")

df$phase <- factor(df$phase, levels=c("rapi", "trig", "shot", "inac"))
df$sampleid <- factor(df$sampleid, levels=c("p30_0.1", "p30_1", "p30_10","p30_100"))

g<-ggplot()+
  geom_line(data=df,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generations")+
  ylab("insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  ggtitle("piRNA clusters = 0%, paramutable loci = 30%") +
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                                      c("p30_0.1" = "trigger loci = 0.1%",
                                                        "p30_1" = "Trigger loci = 1%",
                                                        "p30_10" = "Trigger loci = 10%",
                                                        "p30_100" = "Trigger loci = 100%")))

plot(g)
```

![](2022_08_09_Simulation_2_Trigger_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
df_2<-read.table("2023_01_17_Simulation_2_Trigger_2", fill = TRUE, sep = "\t")
names(df_2)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
               "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
               "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")


df_2$phase <- factor(df_2$phase, levels=c("rapi", "trig", "shot", "inac"))
df_2$sampleid <- factor(df_2$sampleid, levels=c("p10_0.1", "p10_1", "p10_10","p10_100"))


g_2<-ggplot()+
  geom_line(data=df_2,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generations")+
  ylab("insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  ggtitle("piRNA clusters = 0%, paramutable loci = 10%") +
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("p10_0.1" = "Trigger loci = 0.1%",
                                                "p10_1" = "Trigger loci = 1%",
                                                "p10_10" = "Trigger loci = 10%",
                                                "p10_100" = "Trigger loci = 100%")))

plot(g_2)
```

![](2022_08_09_Simulation_2_Trigger_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
df_3<-read.table("2023_01_17_Simulation_2_Trigger_3", fill = TRUE, sep = "\t")
names(df_3)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq",
               "fixed","spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4",
               "fwpar_yespi","fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")


df_3$phase <- factor(df_3$phase, levels=c("rapi", "trig", "shot", "inac"))
df_3$sampleid <- factor(df_3$sampleid, levels=c("p3_10_0.1", "p3_10_1", "p3_10_10","p3_10_100"))


g_3<-ggplot()+
  geom_line(data=df_3,aes(x=gen,y=avtes,group=rep,color=phase),alpha=1,size=0.7)+
  xlab("generations")+
  ylab("insertions per diploid individual")+
  theme(legend.position="none")+
  scale_colour_manual(values=p)+
  ylim(0,500)+
  ggtitle("piRNA clusters = 3%, paramutable loci = 10%") +
  facet_wrap(~sampleid, labeller = labeller(sampleid = 
                                              c("p3_10_0.1" = "Trigger loci = 0.1%",
                                                "p3_10_1" = "Trigger loci = 1%",
                                                "p3_10_10" = "Trigger loci = 10%",
                                                "p3_10_100" = "Trigger loci = 100%")))

plot(g_3)
```

![](2022_08_09_Simulation_2_Trigger_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

## Conclusions

1% of the genome as trigger loci is enough to reduce the average number
of transposable elements both with 10% and 30% paramutable loci.
