Validation piRNA clusters
================
Almorò Scarpa

## Introduction

### Scenarios

-   pc0: no piRNA clusters

-   pc50: 50% of the genome is composed by piRNA clusters

-   pc100: 100% of the genome is composed by piRNA clusters

### Starting conditions equal for all the simulations

![c\_{0}=10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7B0%7D%3D10 "c_{0}=10")

![\\mu = 0.1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%20%3D%200.1 "\mu = 0.1")

![t=100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t%3D100 "t=100")

A chromosome of size 1Mb and no piRNA clusters. 500 replicates.

## Materials & Methods

version: invadego0.2.1

-   seed pc0: 01659107659347760000

-   seed pc50: 1659107659391276000

-   seed pc100: 1659107659369077000

### Commands for the simulation:

``` bash
folder="/Users/ascarpa/Paramutations_TEs/Validation/Raw"
tool="/Users/ascarpa/invade-invadego/invadego021"

$tool --N 1000 --gen 100 --genome mb:1 --cluster kb:0 --rr 4 --rep 100 --u 0.1 --basepop 10 --steps 1 --sampleid pc0> $folder/validation_3_1_clusters &
$tool --N 1000 --gen 100 --genome mb:1 --cluster kb:500 --rr 4 --rep 100 --u 0.1 --basepop 10 --steps 1 --sampleid pc50> $folder/validation_3_2 &
$tool --N 1000 --gen 100 --genome mb:1 --cluster kb:1000 --rr 4 --rep 100 --u 0.1 --basepop 10 --steps 1 --sampleid pc100> $folder/validation_3_3

cat validation_3_1 validation_3_2 validation_3_3 |grep -v "^Invade"|grep -v "^#" > 2022_08_01_Validation_3_piRNA_clusters
```

### Visualization in R

Setting the environment

``` r
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(gridExtra)
theme_set(theme_bw())
```

Visualization: comparing the simulations with the prediction

``` r
validation<-read.table("Raw/2022_08_01_Validation_3_piRNA_clusters", fill = TRUE, sep = "\t")
names(validation)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi","avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid")
data_new <- validation
data_new$sampleid <- factor(data_new$sampleid,
                         levels = c("pc0", "pc50", "pc100"))

gl<-ggplot()+geom_line(data=data_new,aes(x=gen,group=rep,y=avtes*1000),alpha=0.4)+scale_y_log10()+theme(legend.position="none")+ylab("log10 TE copies in population")+xlab("generation")+facet_grid(.~sampleid)
plot(gl)
```

![](2022_08_01_Validation_3_piRNA_clusters_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Conclusions

The simulation matched the expectations.

The introduction of piRNA clusters stops the spreading of TEs in the
population.
