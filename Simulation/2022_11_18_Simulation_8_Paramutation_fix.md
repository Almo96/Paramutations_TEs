2022_11_18_Simulation_8\_Paramutation_fix
================
Almo
2022-11-18

## Introduction

With this simulation we wanted to study the possibility of fixation of a
paramutation in a population.

### Initial conditions:

A population of 1000, 1 chromosomes of size 10 Mb, no piRNA clusters and
a TE equal of for each chromosome in position 1, half of the population
with maternal piRNA deposition.

We used 1000 replicates for each simulation.

## Materials & Methods

version: invadego0.23

-   seed 8_1:

-   seed 8_2:

-   seed 8_3:

-   seed 8_4:

-   seed 8_5:

### Commands for the simulation:

``` bash
folder="/Users/ascarpa/Paramutations_TEs/Simulation/Raw"
tool="/Users/ascarpa/invade-invadego/invadego023"

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0 --gen 5000 --genome mb:10 --steps 5000 --rr 4 --paramutation 999999:1 --rep 1000 --silent > $folder/2022_11_18_simulation_8_1

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.01 -x 0.01 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_2

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.01 -x 0.1 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_3

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.1 -x 0.01 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_4

$tool --N 1000 --basepop $folder/2022_11_18_input_08 --cluster kb:0 --u 0.1 -x 0.1 --gen 500 --genome mb:10 --steps 100 --rr 4 --paramutation 999999:1 --rep 100 --silent > $folder/2022_11_18_simulation_8_5
```

### Visualization in R

Setting the environment

Visualization:

## Conclusions
