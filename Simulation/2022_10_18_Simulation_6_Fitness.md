2022_10_18_Simulation_6\_Fitness
================
Almo
2022-10-18

## Introduction

With this simulation we wanted to understand the role of paramutations
on fitness during a TEs invasion.

### Initial conditions:

A population of 1000, 5 chromosomes of size 10 Mb, 5 piRNA clusters of
size 300 Kb and an initial number of TEs in the population equal to 100.

We used 100 replicates for each simulation.

## Materials & Methods

version: invadego0.2.2

### Commands for the simulation:

``` bash
folder="/Users/ascarpa/Paramutations_TEs/Simulation/Raw"
tool="/Users/ascarpa/invade-invadego/invadego022"

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --steps 1 --sampleid p0_x0 > $folder/2022_10_18_simulation_6_1 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 --steps 1 --sampleid p0_x0.02 > $folder/2022_10_18_simulation_6_2 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 --steps 1 --sampleid p0_x0.03 > $folder/2022_10_18_simulation_6_3 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.05 --steps 1 --sampleid p0_x0.05 > $folder/2022_10_18_simulation_6_4 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 100:1 --steps 1 --sampleid p1_x0 > $folder/2022_10_18_simulation_6_5 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 --paramutation 100:1 --steps 1 --sampleid p1_x0.02 > $folder/2022_10_18_simulation_6_6 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 --paramutation 100:1 --steps 1 --sampleid p1_x0.03 > $folder/2022_10_18_simulation_6_7 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.5 --basepop 100 -x 0.05 --paramutation 100:1 --steps 1 --sampleid p1_x0.05 > $folder/2022_10_18_simulation_6_8 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 --steps 1 --sampleid p10_x0 > $folder/2022_10_18_simulation_6_9 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 --paramutation 10:1 --steps 1 --sampleid p10_x0.02 > $folder/2022_10_18_simulation_6_10 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 --paramutation 10:1 --steps 1 --sampleid p10_x0.03 > $folder/2022_10_18_simulation_6_11 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.05 --paramutation 10:1 --steps 1 --sampleid p10_x0.05 > $folder/2022_10_18_simulation_6_12 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 1:0 --steps 1 --sampleid p100_x0 > $folder/2022_10_18_simulation_6_13 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 --paramutation 1:0 --steps 1 --sampleid p100_x0.02 > $folder/2022_10_18_simulation_6_14 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 --paramutation 1:0 --steps 1 --sampleid p100_x0.03 > $folder/2022_10_18_simulation_6_15 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.05 --paramutation 1:0 --steps 1 --sampleid p100_x0.05 > $folder/2022_10_18_simulation_6_16 &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -no-x-cluins --steps 1 --sampleid p0_x0 > $folder/2022_10_18_simulation_6_1_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 -no-x-cluins --steps 1 --sampleid p0_x0.02 > $folder/2022_10_18_simulation_6_2_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 -no-x-cluins --steps 1 --sampleid p0_x0.03 > $folder/2022_10_18_simulation_6_3_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.05 -no-x-cluins --steps 1 --sampleid p0_x0.05 > $folder/2022_10_18_simulation_6_4_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 100:1 -no-x-cluins --steps 1 --sampleid p1_x0 > $folder/2022_10_18_simulation_6_5_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 --paramutation 100:1 -no-x-cluins --steps 1 --sampleid p1_x0.02 > $folder/2022_10_18_simulation_6_6_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 --paramutation 100:1 -no-x-cluins --steps 1 --sampleid p1_x0.03 > $folder/2022_10_18_simulation_6_7_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.5 --basepop 100 -x 0.05 --paramutation 100:1 -no-x-cluins --steps 1 --sampleid p1_x0.05 > $folder/2022_10_18_simulation_6_8_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 10:1 -no-x-cluins --steps 1 --sampleid p10_x0 > $folder/2022_10_18_simulation_6_9_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 --paramutation 10:1 -no-x-cluins --steps 1 --sampleid p10_x0.02 > $folder/2022_10_18_simulation_6_10_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 --paramutation 10:1 -no-x-cluins --steps 1 --sampleid p10_x0.03 > $folder/2022_10_18_simulation_6_11_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.05 --paramutation 10:1 -no-x-cluins --steps 1 --sampleid p10_x0.05 > $folder/2022_10_18_simulation_6_12_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 --paramutation 1:0 -no-x-cluins --steps 1 --sampleid p100_x0 > $folder/2022_10_18_simulation_6_13_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.02 --paramutation 1:0 -no-x-cluins --steps 1 --sampleid p100_x0.02 > $folder/2022_10_18_simulation_6_14_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.03 --paramutation 1:0 -no-x-cluins --steps 1 --sampleid p100_x0.03 > $folder/2022_10_18_simulation_6_15_noxclu &

$tool --N 1000 --gen 5000 --genome mb:10,10,10,10,10 --cluster kb:300,300,300,300,300 --rr 4,4,4,4,4 --rep 100 --u 0.1 --basepop 100 -x 0.05 --paramutation 1:0 -no-x-cluins --steps 1 --sampleid p100_x0.05 > $folder/2022_10_18_simulation_6_16_noxclu

cat 2022_10_18_simulation_6_1 2022_10_18_simulation_6_5 2022_10_18_simulation_6_9 2022_10_18_simulation_6_13 |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0

cat 2022_10_18_simulation_6_2 2022_10_18_simulation_6_6 2022_10_18_simulation_6_10  2022_10_18_simulation_6_14 |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0.02

cat 2022_10_18_simulation_6_3 2022_10_18_simulation_6_7 2022_10_18_simulation_6_11  2022_10_18_simulation_6_15 |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0.03

cat 2022_10_18_simulation_6_4 2022_10_18_simulation_6_8 2022_10_18_simulation_6_12  2022_10_18_simulation_6_16 |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0.05

cat 2022_10_18_simulation_6_1_noxclu 2022_10_18_simulation_6_5_noxclu 2022_10_18_simulation_6_9_noxclu 2022_10_18_simulation_6_13_noxclu |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0_noxclu

cat 2022_10_18_simulation_6_2_noxclu 2022_10_18_simulation_6_6_noxclu 2022_10_18_simulation_6_10_noxclu  2022_10_18_simulation_6_14_noxclu |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0.02_noxclu

cat 2022_10_18_simulation_6_3_noxclu 2022_10_18_simulation_6_7_noxclu 2022_10_18_simulation_6_11_noxclu  2022_10_18_simulation_6_15_noxclu |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0.03_noxclu

cat 2022_10_18_simulation_6_4_noxclu 2022_10_18_simulation_6_8_noxclu 2022_10_18_simulation_6_12_noxclu  2022_10_18_simulation_6_16_noxclu |grep -v "^Invade"|grep -v "^#" > 2022_10_18_Simulation_6_Fitness_0.05_noxclu
```

### Visualization in R

Setting the environment

``` r
library(tidyverse)
library(ggplot2)
theme_set(theme_bw())
```

Visualization:

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-5.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-6.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-7.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-8.png)<!-- -->

``` r
g_A_0.03_noxclu <- ggplot(df_0.03_noxclu,aes(x = gen, y = avtes))+
  geom_line(alpha=0.1)+
  xlab("generation")+
  ylab("TEs insertions per diploid individual")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid = 
                                              c("p0_x0.03" = "Paramutable loci = 0%",
                                                "p1_x0.03" = "Paramutable loci = 1%",
                                                "p10_x0.03" = "Paramutable loci = 10%",
                                                "p100_x0.03" = "Paramutable loci = 100%")))
plot(g_A_0.03_noxclu)
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-9.png)<!-- -->

``` r
g_B_0.03_noxclu <- ggplot(df_0.03_noxclu,aes(x = gen, y = avw))+
  geom_line(alpha=0.12)+
  xlim(0,2500)+
  xlab("generation")+
  ylab("Fitness")+
  facet_wrap(~sampleid, ncol=4, labeller = labeller(sampleid =
                                                      c("p0_x0.03" = "Paramutable loci = 0%",
                                                        "p1_x0.03" = "Paramutable loci = 1%",
                                                        "p10_x0.03" = "Paramutable loci = 10%",
                                                        "p100_x0.03" = "Paramutable loci = 100%")))
plot(g_B_0.03_noxclu)
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-10.png)<!-- -->

``` r
df_summary <- df_0.03_noxclu %>% 
  group_by(sampleid, rep) %>%
  summarize(min_fitness = min(avw))

df_summary$sampleid[df_summary$sampleid == "p0_x0.03"] <- "0% (Trap model)"
df_summary$sampleid[df_summary$sampleid == "p1_x0.03"] <- "1%"
df_summary$sampleid[df_summary$sampleid == "p10_x0.03"] <- "10%"
df_summary$sampleid[df_summary$sampleid == "p100_x0.03"] <- "100%"

boxplot(1-(df_summary$min_fitness) ~ df_summary$sampleid, xlab="Paramutable loci", ylab = "Fitness cost")
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-11.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-12.png)<!-- -->

``` r
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
```

![](2022_10_18_Simulation_6_Fitness_files/figure-gfm/unnamed-chunk-3-13.png)<!-- -->

## Conclusions

An increase in paramutations reduces the fitness cost of the TEs during
an invasion.
