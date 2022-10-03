Validation of invasion
================
Almorò Scarpa

## Introduction

With this validation we wanted to prove that the simulator will generate
on average the number of insertions predicted by the equation:

![c\_{t}=c\_{0}(1+\\mu)^t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7Bt%7D%3Dc_%7B0%7D%281%2B%5Cmu%29%5Et "c_{t}=c_{0}(1+\mu)^t")

![c\_{t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7Bt%7D "c_{t}")
TE copies at generation t

![c\_{0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7B0%7D "c_{0}")
TE copies at generation 0

![\\mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu "\mu")
Transposition rate

![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
Generation

### Initial conditions:

![c\_{0}=10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;c_%7B0%7D%3D10 "c_{0}=10")

![\\mu = 0.1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmu%20%3D%200.1 "\mu = 0.1")

![t=100](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t%3D100 "t=100")

A chromosome of size 1Mb and no piRNA clusters. We used 500 replicates

## Materials & Methods

version: invadego0.2.1 seed: 1659093816078434000

### Commands for the simulation:

``` bash
folder="/Users/ascarpa/Paramutations_TEs/Validation/Raw"
tool="/Users/ascarpa/invade-invadego/invadego021"

$tool --N 1000 --gen 100 --genome mb:1 --cluster kb:0 --rr 4 --rep 500 --u 0.1 --basepop 10 --silent --steps 1 > $folder/2022_08_01_Validation_1_invasion
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
cn<-seq(0,99,1)
res<-10*1.1^cn
theo<-data.frame(x=1:100,y=res)
validation<-read.table("Raw/2022_08_01_Validation_1_invasion")
names(validation)<-c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "avtes", "avpopfreq", "fixed","spacer_2","phase","fwpirna","spacer_3","fwcli","avcli","fixcli","spacer_4","fwpar_yespi","fwpar_nopi",
                     "avpar","fixpar","spacer_5","piori")

gl<-ggplot()+geom_line(data=validation,aes(x=gen,group=rep,y=avtes*1000),alpha=0.15,size=0.3)+scale_y_log10()+geom_line(data=theo,aes(x=x,y=y),size=2)+theme(legend.position="none")+ylab("log10 TE copies in population")+xlab("Generation")
plot(gl)
```

![](2022_08_01_Validation_1_invasion_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Conclusions

The simulation matched the expectations. Invadego accurately reproduces
the expected exponential increase of TE copy numbers in a population
with no piRNA clusters.
