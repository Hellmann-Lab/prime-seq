
## Purpose: Figure 2 – PCA

### 1\. Load the following packages:

``` r
library(tidyverse)
library(ggsignif)
library(ggrepel)
library(edgeR)
library(genefilter)
library(grid)
library(gridExtra)
library(ggsci)
library(dplyr)
library(UpSetR)
library(cowplot)
library(Seurat)
```

### 2\. Load following functions:

``` r
## all necessary custom functions are in the following script
source(paste0(here::here(),"/0_Scripts/custom_functions.R"))

theme_pub <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"),
                                     axis.text = element_text(colour="black", size=14), 
                                     axis.title=element_text(size=16,face="bold"), 
                                     legend.text=element_text(size=14),
                                     legend.position="right",
                                     axis.line.x = element_line(colour = "black"), 
                                     axis.line.y = element_line(colour = "black"),
                                     strip.background=element_blank(), 
                                     strip.text=element_text(size=16))  

#prevent scientific notation
options(scipen=999)
fig_path<-paste0(here::here(),"/1_RNA_isolation/")
```

### 3\. HEK - Load Data

``` r
counts <- readRDS(paste0(fig_path,"Bulk_opt_lysis_test_2_HEK.dgecounts.rds"))

inf <- read.csv(paste0(fig_path,"sample_info.csv"), header = T, stringsAsFactors = F)

inf$Sample <- as.character(inf$Sample)

inf<-inf %>% 
  mutate(Condition=case_when(Condition=="Incubation + ProtK"~"Magnetic Beads",
                             TRUE~Condition))

inf_HEK <- inf %>% filter(Celltype == "HEK") %>% filter(Cells == "10")

#inex
inex_umi <- as.matrix(counts$umicount$inex$all) %>% remove_Geneversion()
```

### 4\. HEK - Clustering

``` r
inex_umi <- inex_umi[,inf_HEK$BC]
inex_umi <- inex_umi[grep(rownames(inex_umi),pattern="ERCC*",invert = T),]
data_seurat <- CreateSeuratObject(counts = inex_umi, project = "HEK", min.cells = 3, min.features = 200)
data_seurat <- NormalizeData(data_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data_seurat)
data_seurat <- ScaleData(data_seurat, features = all.genes)
data_seurat@meta.data$condition = inf_HEK$Condition

mat<-t(data_seurat@assays$RNA@scale.data)
pc<-prcomp(mat, scale=F)
#inf_HEK = inf_HEK[rownames(pc$x), ]
percVar <- data.frame(pc$sdev^2/sum(pc$sdev^2)*100)
pcs<-data.frame(pc$x,condition = factor(inf_HEK$Condition, levels = c("Magnetic Beads", "Incubation","No Incubation","Column")))

PCA12 <- ggplot(data= pcs, aes(x=PC1, y=PC2, col = condition)) +
  geom_point(size =4, aes(shape=condition))+
  theme_pub+
  ylim(-100,100)+
  xlim(-150, 150)+
  scale_color_manual(values = c("#008080","gray70", "gray60", "gray40"))+
  theme(legend.position="bottom",
        legend.title = element_blank())+
  xlab(paste0("PC1: ",round(percVar[1,1]),"% variance")) +
  ylab(paste0("PC2: ",round(percVar[2,1]),"% variance"))

PCA12
```

![](Lysis_PCA_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
PCA23 <- ggplot(data= pcs, aes(x=PC2, y=PC3, col = condition)) +
  geom_point(size =4, aes(shape=condition))+
  theme_pub+
  scale_color_manual(values = c("#008080","gray70", "gray60", "gray40"))+
  theme(legend.position="bottom")+
  xlab(paste0("PC2: ",round(percVar[2,1]),"% variance")) +
  ylab(paste0("PC3: ",round(percVar[3,1]),"% variance"))

PCA23
```

![](Lysis_PCA_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->