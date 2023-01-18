library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(patchwork)

setwd("/Users/ascarpa/Paramutations_TEs/Simulation/Raw")
df<-read.table("2022_11_07_Simulation_7_sel_para_u02_clu100", sep = "\t", fill = TRUE, row.names=NULL)
df_2<-read.table("2022_11_7_Simulation_7_clu_para_u02_x001", sep = "\t", fill = TRUE, row.names=NULL)
df_para_0_sel<-read.table("2022_11_30_Simulation_7_para0_sel", sep = "\t", fill = TRUE, row.names=NULL)
df_para_0_clu<-read.table("2022_11_30_Simulation_7_para0_clu", sep = "\t", fill = TRUE, row.names=NULL)

naming <- c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq", "fixed",
            "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
            "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid_x", "sampleid_para","extra")

naming_2 <- c("rep", "gen", "popstat", "fmale", "spacer_1", "fwte", "avw", "min_w", "avtes", "avpopfreq", "fixed",
              "spacer_2", "phase", "fwpirna", "spacer_3", "fwcli", "avcli", "fixcli", "spacer_4", "fwpar_yespi",
              "fwpar_nopi", "avpar","fixpar","spacer_5","piori","orifreq","spacer 6", "sampleid_clu", "sampleid_para","extra")


names(df) <- naming
names(df_2) <- naming_2
names(df_para_0_sel) <- naming[-c(30)]
names(df_para_0_clu) <- naming_2[-c(30)]


r <- runif(1000, -1, 1)
df_para_0_sel$sampleid_para <- r
df_para_0_clu$sampleid_para <- r


df$sampleid_x <- as.numeric(df$sampleid_x)
df$sampleid_para <- as.numeric(df$sampleid_para)
df$min_w <- as.numeric(df$min_w)

df_2$sampleid_clu <- as.numeric(df_2$sampleid_clu)
df_2$sampleid_para <- as.numeric(df_2$sampleid_para)
df_2$min_w <- as.numeric(df_2$min_w)

df_para_0_sel$sampleid_x <- as.numeric(df_para_0_sel$sampleid_x)
df_para_0_sel$sampleid_para <- as.numeric(df_para_0_sel$sampleid_para)
df_para_0_sel$min_w <- as.numeric(df_para_0_sel$min_w)

df_para_0_clu$sampleid_x <- as.numeric(df_para_0_clu$sampleid_clu)
df_para_0_clu$sampleid_para <- as.numeric(df_para_0_clu$sampleid_para)
df_para_0_clu$min_w <- as.numeric(df_para_0_clu$min_w)

#Keep only last generation, will be less then 5000 if fail
df<-subset(df, gen > 0)
df_2<-subset(df_2, gen > 0)
df_para_0_sel<-subset(df_para_0_sel, gen > 0)
df_para_0_clu<-subset(df_para_0_clu, gen > 0)

median_w=(median(df_2$min_w))
mean_w=(max(df_2$min_w)+min(df_2$min_w))/2
df_3 <- subset(df_2, select=c("rep","min_w","sampleid_para"))
df_3$d_m <- sqrt((df_3$min_w-median_w)**2)

d_points <- df_3 %>%
  dplyr::group_by(sampleid_para) %>%
  dplyr::summarize(min_d_m=min(d_m))

df_3 <- na.omit(df_3)
d_points <- na.omit(d_points)

df_3 <-df_3[order(df_3$sampleid_para),]

df_4 <- as.data.frame(matrix(nrow=500,ncol=3))
colnames(df_4)<-c("para","min_w", "rep")



x=0
x_1=0
y=1
k=1
while (x<nrow(df_4)-1){
  x_1=x+1
  df_4[x_1,1]=x
  y=k
  while(y<nrow(df_3)){
    if (d_points[x_1,1] == df_3[y,3] && d_points[x_1,2] == df_3[y,4]){
      df_4[x_1,2]<-df_3[y,2]
      df_4[x_1,3]<-df_3[y,1]
    }
    if (d_points[x_1,1] < df_3[y,3]){
      k=y
      y=nrow(df_3)
    }
    y=y+1
  }
  x=x+1
}


rep_num <- as.vector(df_4$rep)
rep_num <- sort(rep_num)
df_2_cleaned <- df_2[df_2$rep %in% rep_num, ]


#Unique gradient
min_df=min(df$min_w)
min_sel_para0=min(df_para_0_sel$min_w)

if(min_df < min_sel_para0){
  min_sel = min_df
} else {
  min_sel = min_sel_para0
}


min_df2=min(df_2$min_w)
min_clu_para0=min(df_para_0_clu$min_w)

if(min_df2 < min_clu_para0){
  min_clu = min_df2
} else {
  min_clu = min_clu_para0
}


color.gradient_sel <- function(x, colors=c("#D7191C","#FDAE61","#A6D96A","#1A9641"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min_sel,1.0, length.out=colsteps)) ] )
}

color.gradient_clu <- function(x, colors=c("#D7191C","#FDAE61","#A6D96A","#1A9641"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min_clu,1.0, length.out=colsteps)) ] )
}

df$col<-color.gradient_sel(df$min_w)
df[df$popstat=="fail-w",]$col<-"white"
df[df$popstat=="fail-0",]$col<-"grey"
df$col<-as.factor(df$col)

df_2$col<-color.gradient_clu(df_2$min_w)
df_2[df_2$popstat=="fail-w",]$col<-"white"
df_2[df_2$popstat=="fail-0",]$col<-"grey"
df_2$col<-as.factor(df_2$col)


df_para_0_sel$col<-color.gradient_sel(df_para_0_sel$min_w)
df_para_0_sel[df_para_0_sel$popstat=="fail-w",]$col<-"white"
df_para_0_sel[df_para_0_sel$popstat=="fail-0",]$col<-"grey"
df_para_0_sel$col<-as.factor(df_para_0_sel$col)

df_para_0_clu$col<-color.gradient_clu(df_para_0_clu$min_w)
df_para_0_clu[df_para_0_clu$popstat=="fail-w",]$col<-"white"
df_para_0_clu[df_para_0_clu$popstat=="fail-0",]$col<-"grey"
df_para_0_clu$col<-as.factor(df_para_0_clu$col)


g_para_selection <- ggplot(df,aes(x=sampleid_para/10,y=sampleid_x,color=col))+scale_color_manual(values=levels(df$col))+
  geom_point(alpha=0.7,size=0.8)+scale_y_log10()+
  ylab("Negative selection coefficient")+
  xlab("Size of paramutable loci as [%]")+
  geom_hline(aes(yintercept=0.001), linetype = "dashed", size=1)+
  theme(legend.position = "none",panel.background = element_rect(fill="grey90"))

plot(g_para_selection)


pdf(file = "/Users/ascarpa/Paramutations_TEs/Pictures_paper/Figure_4/para_sel.pdf", width = 10, height = 7.5)
g_para_selection
dev.off()

g_para_cluster<-ggplot(df_2,aes(x=sampleid_para/1000,y=sampleid_clu/100000,color=col))+scale_color_manual(values=levels(df_2$col))+
  geom_point(alpha=0.7,size=0.8)+scale_y_log10()+
  geom_curve(aes(x = 0, y = 1.6, xend = 0.5, yend = 0.30), data = df_2, angle = 90, curvature = 0.2, color="#4B5320", alpha=0.5)+
  scale_x_continuous(labels = scales::percent)+
  ylab("Size of piRNA clusters as [%]")+
  xlab("Size of paramutable loci as [%]")+
  theme(legend.position = "none",panel.background = element_rect(fill="grey90"))

plot(g_para_cluster)


g_para_cluster2 <-ggplot(df_2)+scale_color_manual(values=levels(df_2$col))+
  geom_point(aes(x=sampleid_para/10,y=sampleid_clu/100000,color=col), alpha=0.7,size=0.8)+scale_y_log10()+
  geom_smooth(data=df_2_cleaned,aes(x=sampleid_para/10,y=sampleid_clu/100000),  se = FALSE, method = "loess", formula = y ~ x, colour="black")+
  ylab("Size of piRNA clusters as [%]")+
  xlab("Size of paramutable loci as [%]")+
  theme(legend.position = "none",panel.background = element_rect(fill="grey90"))

plot(g_para_cluster2)

g_para_selection + g_para_cluster2

pdf(file = "/Users/ascarpa/Paramutations_TEs/Pictures_paper/Figure_4/2_graphics_lines.pdf", width = 10, height = 7.5)
g_para_selection + g_para_cluster2
dev.off()


g_para_0_selection <- ggplot(df_para_0_sel,aes(x=sampleid_para/10,y=sampleid_x,color=col))+scale_color_manual(values=levels(df_para_0_sel$col))+
  geom_point(alpha=0.7,size=0.8)+scale_y_log10()+
  ylab("Negative selection coefficient")+
  xlab("Size of paramutable loci as [%]")+
  geom_hline(aes(yintercept=0.001), linetype = "dashed", size=1)+
  theme(legend.position = "none",panel.background = element_rect(fill="grey90"))

plot(g_para_0_selection)


g_para_0_cluster <-ggplot(df_para_0_clu)+scale_color_manual(values=levels(df_para_0_clu$col))+
  geom_point(aes(x=sampleid_para/10,y=sampleid_clu/100000,color=col), alpha=0.7,size=0.8)+scale_y_log10()+
  ylab("Size of piRNA clusters as [%]")+
  xlab("Size of paramutable loci as [%]")+
  theme(legend.position = "none",panel.background = element_rect(fill="grey90"))

plot(g_para_0_cluster)

(g_para_0_selection + g_para_selection) / (g_para_0_cluster + g_para_cluster2)

pdf(file = "/Users/ascarpa/Paramutations_TEs/Pictures_paper/Figure_4/4_graphics.pdf", width = 10, height = 7.5)
(g_para_0_selection + g_para_selection) / (g_para_0_cluster + g_para_cluster2)
dev.off()

(g_para_0_selection | g_para_0_selection | g_para_0_selection | g_para_0_selection | g_para_0_selection | g_para_0_selection) /
  (g_para_0_cluster | g_para_0_cluster | g_para_0_cluster | g_para_0_cluster | g_para_0_cluster | g_para_0_cluster)

pdf(file = "/Users/ascarpa/Paramutations_TEs/Pictures_paper/Figure_4/para_0.pdf", width = 10, height = 7.5)
(g_para_0_selection | g_para_0_selection | g_para_0_selection | g_para_0_selection) /
  (g_para_0_cluster | g_para_0_cluster | g_para_0_cluster | g_para_0_cluster)
dev.off()


color_temp="#EFB462"
column_a_b <- c(0,1,2,3)
column_c <- c(color_temp, color_temp, color_temp, color_temp)

dfab <- data.frame(column_a_b, column_a_b, column_c)


ab <- ggplot(dfab,aes(x = column_a_b, y = column_a_b))+scale_color_manual(values=levels(dfab$column_c))+
  geom_point(alpha=0.7,size=4, color=column_c)

plot(ab)

#D7191C","#FDAE61","#A6D96A","#1A9641"
#> min(df$min_w)
#[1] 0.43
#> max(df$min_w)
#[1] 1
#> color.gradient_sel(0.43)
#[1] "#DD3428"
#> color.gradient_sel(1)
#[1] "#1A9641"
#> (1-0.43)/2
#[1] 0.285
#> 0.285+0.43
#[1] 0.715
#> color.gradient_sel(0.715)
#[1] "#CAC666"


#> min(df_2$min_w)
#[1] 0.1
#> max(df_2$min_w)
#[1] 0.96
#  > color.gradient_clu(0.96)
#[1] "#2FA047"
#  > color.gradient_clu(0.1)
#[1] "#DA2622"

#color.gradient_clu(0.43)
#[1] "#EFB462"
