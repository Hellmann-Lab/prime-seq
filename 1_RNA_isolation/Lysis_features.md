
## Purpose:

Figure 2 - feature plots

## Protocol:

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
library(cowplot)
```

### 2\. Load following functions:

``` r
## all necessary custom functions are in the following script
source("/data/share/htp/prime-seq_Paper/Scripts/custom_functions.R")

theme_pub <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=18, face="bold"),
                                     axis.text = element_text(colour="black", size=14), 
                                     axis.title=element_text(size=16,face="bold"), 
                                     legend.text=element_text(size=14),
                                     legend.position="right",
                                     axis.line.x = element_line(colour = "black"), 
                                     axis.line.y = element_line(colour = "black"),
                                     strip.background=element_blank(), 
                                     strip.text=element_text(size=16))  

theme_set(theme_pub)

#prevent scientific notation
options(scipen=999)

fig_path<-"/data/share/htp/prime-seq_Paper/Fig_beads_columns/"
```

# Features Plots

### 3\. Information files

Read the information table, which also contains the total sequenced
reads calculated using:

zcat file.fq.gz | wc -l

then divide by 4

``` r
feat_seq_info <- read_delim(paste0(fig_path,"analysis/Lysis/Features_sequencing_info.csv"),
                            ";", 
                            escape_double = FALSE, 
                            trim_ws = TRUE)

inf <- read.csv(paste0(fig_path,"analysis/Lysis/sample_info.csv"), header = T, stringsAsFactors = F)

inf$Sample <- as.character(inf$Sample)

inf<-inf %>% 
  mutate(Condition=case_when(Condition=="Incubation + ProtK"~"Magnetic Beads",
                             TRUE~Condition))
```

### 4\. readspercell files

Load readspercell files from zUMIs, which show all features for every i7

``` r
#read files
readspercell_HEK <- read.table(paste0(fig_path,"zUMIs/HEK/zUMIs_output/stats/Bulk_opt_lysis_test_2_HEK.readspercell.txt"), header = T)
readspercell_HEK$Celltype <- "HEK"

readspercell_PBMCs <- read.table(paste0(fig_path,"zUMIs/PBMC/zUMIs_output/stats/Bulk_opt_lysis_PBMCs.readspercell.txt"), header = T)
readspercell_PBMCs$Celltype <- "PBMCs"

readspercell_Tissue <- read.table(paste0(fig_path,"zUMIs/Tissue/zUMIs_output/stats/Bulk_opt_lysis_Tissue.readspercell.txt"), header = T)
readspercell_Tissue$Celltype <- "Tissue"


readspercell <- bind_rows(readspercell_HEK, 
                          readspercell_PBMCs,
                          readspercell_Tissue)

colnames(readspercell)[1] <- "XC" 
```

### 5\. Sequencing Features Bar Plot

Calculate the necessary features for the sequencing table

``` r
#calculate the reads that passed phred QC from zUMIs 
QCPass <- readspercell %>% group_by(Celltype) %>% summarize(QCPass=sum(N)) %>% as.data.frame()

feat_seq_info <- left_join(feat_seq_info, QCPass, by = "Celltype")

#calculate the reads that failed phred QC from zUMIs
feat_seq_info$QCFail <- feat_seq_info$Total_after_filtering - feat_seq_info$QCPass

#calculate unmapped reads that passed QC
Unmapped <- readspercell %>% group_by(Celltype) %>% filter(type == "Unmapped") %>% summarize(Unmapped=sum(N)) %>% as.data.frame()

feat_seq_info <- left_join(feat_seq_info, Unmapped, by = "Celltype")

  
#calculate mapped reads from unassigned bcs
Unassigned <- readspercell %>% group_by(Celltype) %>% filter(type != "Unmapped") %>% filter(XC == "bad") %>% summarize(Unassigned=sum(N)) %>% as.data.frame()

feat_seq_info <- left_join(feat_seq_info, Unassigned, by = "Celltype")

#calculate the mapped reads that are assigned to BCs
feat_seq_info$Assigned <- feat_seq_info$QCPass - feat_seq_info$Unmapped - feat_seq_info$Unassigned

#make into long table
feat_seq_info_long <- tidyr::gather(feat_seq_info, type, measurement, QCFail:Assigned, factor_key=TRUE)

#change levels of data
feat_seq_info_long$type<-factor(feat_seq_info_long$type,levels=rev(c("Assigned", "Unassigned", "Unmapped","QCFail")))


#plot sequencing features
plot_feat_seq <- ggplot(subset(feat_seq_info_long, type != "QCFail"), aes(x=1, y=measurement, fill=type))+
    geom_bar(stat="identity", position = "fill")+
    facet_wrap(~Celltype, ncol = 1, strip.position = c("right")) +
    ylab("Fraction of Seq Reads") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    coord_flip() +
    scale_fill_manual(values = c("gray33", "orangered1","dodgerblue4"))+
    guides(fill = guide_legend(reverse=T))+
    theme(legend.title = element_blank(), 
          legend.position="bottom", 
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(),
          strip.text.y = element_text(angle = 360),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank()) 
plot_feat_seq
```

![](Lysis_features_files/figure-gfm/feat_seq_calc-1.png)<!-- -->

### 6\. Mitochondrial, Ribosomal, and lncRNA Reads - HEK - 10k cells

Calculate the mitochondrial and ribosomal reads from exonic and intronic
count matrix and use this for features plots

``` r

## define colours
iso_type_cols<-c("#00BECF","#EB343C","#F49D5B", "#4BAF66")
names(iso_type_cols)<-c("Column", 
                       "No Incubation",
                       "Incubation",
                       "Magnetic Beads")

seqType_cols<-c("#00BECF","#EB343C","#F49D5B", "#4BAF66")
names(seqType_cols)<-c("Mitochondrial", 
                       "Ribosomal",
                       "lncRNA",
                       "Coding")


#load count matrix
counts <- readRDS(paste0(fig_path,"zUMIs/HEK/zUMIs_output/expression/Bulk_opt_lysis_test_2_HEK.dgecounts.rds"))

#subset exonic and intronic reads
exon <- as.matrix(counts$readcount$exon$all) %>% remove_Geneversion()
intron <- as.matrix(counts$readcount$intron$all) %>% remove_Geneversion()
inex <- as.matrix(counts$readcount$inex$all) %>% remove_Geneversion()


gtype <- data.frame(species = "human", getbiotype("hsapiens_gene_ensembl", species = "human"))

## Summarize exonic counts by biotype and condition

exon_gtype_sum<-exon %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="BC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "No Incubation",
                                               "Incubation",
                                               "Magnetic Beads")))

## Summarize intronic counts by biotype and condition

intron_gtype_sum<-intron %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="BC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "No Incubation",
                                               "Incubation",
                                               "Magnetic Beads")))

## Summarize intronic counts by biotype and condition

inex_gtype_sum<-inex %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="BC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "No Incubation",
                                               "Incubation",
                                               "Magnetic Beads")))

plot_gtype_bar_HEK_all_inex <- ggplot(inex_gtype_sum, aes(x=Condition, y=Counts))+
    geom_col(aes(fill=type),position=position_fill()) +
    ylab("Counts") +
    xlab("Extraction Conditions") +
    scale_fill_manual(values = seqType_cols)+
    facet_grid(Cells~.)+
    coord_flip()+
    theme(strip.background =element_rect(fill="gray90"),
          legend.title = element_blank(),
          legend.position= "right",
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=45,hjust=1))
plot_gtype_bar_HEK_all_inex
```

![](Lysis_features_files/figure-gfm/mito_ribo%20HEK-1.png)<!-- -->

``` r

exon_gtype_sum$seq_type<-"Exon"
intron_gtype_sum$seq_type<-"Intron"
inex_gtype_sum$seq_type<-"Intron Exon"

gtype_sum_all<-bind_rows(exon_gtype_sum,intron_gtype_sum,inex_gtype_sum)

mitoribo_agg_HEK_ex <- exon_gtype_sum %>% 
  group_by(type,cond) %>% 
  summarize(Counts=sum(Counts),Celltype=unique(Celltype),Condition=unique(Condition)) 

mitoribo_agg_HEK_inex <- inex_gtype_sum %>% 
  group_by(type,cond) %>% 
  summarize(Counts=sum(Counts),Celltype=unique(Celltype),Condition=unique(Condition)) 
```

### 7\. Mitochondrial, Ribosomal, and lncRNA Reads - PBMCs

Calculate the mitochondrial and ribosomal reads from exonic and intronic
count matrix and use this for features plots

``` r
#load count matrix
counts <- readRDS(paste0(fig_path,"zUMIs/PBMC/zUMIs_output/expression/Bulk_opt_lysis_PBMCs.dgecounts.rds"))

#subset exonic and intronic reads
exon <- as.matrix(counts$readcount$exon$all) %>% remove_Geneversion()
intron <- as.matrix(counts$readcount$intron$all) %>% remove_Geneversion()
inex <- as.matrix(counts$readcount$inex$all) %>% remove_Geneversion()


gtype <- data.frame( species="human", getbiotype("hsapiens_gene_ensembl", species = "human"))

## Summarize exonic counts by biotype and condition

exon_gtype_sum<-exon %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="XC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column",
                                               "Magnetic Beads")))

## Summarize intronic counts by biotype and condition

intron_gtype_sum<-intron %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="XC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "Magnetic Beads")))

## Summarize intronic counts by biotype and condition

inex_gtype_sum<-inex %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="XC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "Magnetic Beads")))

plot_gtype_bar_PBMC_all_inex <- ggplot(inex_gtype_sum, aes(x=Condition, y=Counts))+
    geom_col(aes(fill=type),position=position_fill()) +
    ylab("Counts") +
    xlab("Extraction Conditions") +
    scale_fill_manual(values = seqType_cols)+
    facet_grid(Cells~.)+
    coord_flip()+
    theme(strip.background =element_rect(fill="gray90"),
          legend.title = element_blank(),
          legend.position= "right",
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=45,hjust=1))
plot_gtype_bar_PBMC_all_inex
```

![](Lysis_features_files/figure-gfm/mito_ribo%20PBMCs-1.png)<!-- -->

``` r

exon_gtype_sum$seq_type<-"Exon"
intron_gtype_sum$seq_type<-"Intron"
inex_gtype_sum$seq_type<-"Intron Exon"

gtype_sum_all<-bind_rows(exon_gtype_sum,intron_gtype_sum,inex_gtype_sum)

mitoribo_agg_PBMC_ex <- exon_gtype_sum %>% 
  group_by(type,cond) %>% 
  summarize(Counts=sum(Counts),Celltype=unique(Celltype),Condition=unique(Condition)) 

mitoribo_agg_PBMC_inex <- inex_gtype_sum %>% 
  group_by(type,cond) %>% 
  summarize(Counts=sum(Counts),Celltype=unique(Celltype),Condition=unique(Condition)) 
```

### 8\. Mitochondrial, Ribosomal, and lncRNA Reads - Tissue

Calculate the mitochondrial and ribosomal reads from exonic and intronic
count matrix and use this for features plots

``` r
#load count matrix
counts <- readRDS(paste0(fig_path,"zUMIs/Tissue/zUMIs_output/expression/Bulk_opt_lysis_Tissue.dgecounts.rds"))

#subset exonic and intronic reads
exon <- as.matrix(counts$readcount$exon$all) %>% remove_Geneversion()
intron <- as.matrix(counts$readcount$intron$all) %>% remove_Geneversion()
inex <- as.matrix(counts$readcount$inex$all) %>% remove_Geneversion()

gtype <- data.frame( species="mouse", getbiotype("mmusculus_gene_ensembl", species = "mouse"))

## Summarize exonic counts by biotype and condition

exon_gtype_sum<-exon %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="XC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "Magnetic Beads")))

## Summarize intronic counts by biotype and condition

intron_gtype_sum<-intron %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="XC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "Magnetic Beads")))

## Summarize intronic counts by biotype and condition

inex_gtype_sum<-inex %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ENSEMBL") %>% 
  left_join(gtype) %>% 
  group_by(type) %>% 
  summarize(across(where(is.double),sum)) %>% 
  mutate(type=if_else(is.na(type),"coding",type)) %>% 
  pivot_longer(cols = where(is.double),values_to="Counts",names_to="XC") %>% 
  mutate(type=case_when(type=="mito"~"Mitochondrial",
                        type=="rrna"~"Ribosomal",
                        type=="lnc"~"lncRNA",
                        type=="coding"~"Coding")) %>% 
  left_join(inf) %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  mutate(Condition=factor(Condition,levels = c("Column", 
                                               "Magnetic Beads")))

plot_gtype_bar_Tissue_all_inex <- ggplot(inex_gtype_sum, aes(x=Condition, y=Counts))+
    geom_col(aes(fill=type),position=position_fill()) +
    ylab("Counts") +
    xlab("Extraction Conditions") +
    scale_fill_manual(values = seqType_cols)+
    coord_flip()+
    theme(strip.background =element_rect(fill="gray90"),
          legend.title = element_blank(),
          legend.position= "right",
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=45,hjust=1))
plot_gtype_bar_Tissue_all_inex
```

![](Lysis_features_files/figure-gfm/mito_ribo%20Tissue-1.png)<!-- -->

``` r
## combine all dfs
exon_gtype_sum$seq_type<-"Exon"
intron_gtype_sum$seq_type<-"Intron"
inex_gtype_sum$seq_type<-"Intron Exon"

gtype_sum_all<-bind_rows(exon_gtype_sum,intron_gtype_sum,inex_gtype_sum)

mitoribo_agg_Tissue_ex <- exon_gtype_sum %>% 
  group_by(type,cond) %>% 
  summarize(Counts=sum(Counts),Celltype=unique(Celltype),Condition=unique(Condition)) 

mitoribo_agg_Tissue_inex <- inex_gtype_sum %>% 
  group_by(type,cond) %>% 
  summarize(Counts=sum(Counts),Celltype=unique(Celltype),Condition=unique(Condition))
```

### 9\. Assigned-Mapped Feature Bar Plot

#### 9.1 Make Table

Calculate the necessary features for the mapped table

``` r
#create a table for mapped features
feat_map_info <- feat_seq_info[,c(1, 7)]

#add condition information to reads per cell
inf<-inf %>% 
  mutate(XC=if_else(Celltype=="HEK",BC,XC))

readspercell_sum <- left_join(readspercell, inf[,3:8], by = "XC") %>% 
  mutate(cond=paste(Celltype, Condition)) %>% 
  filter(XC!="bad")

#aggregate mapped assigned features
mapped <- readspercell_sum %>% 
  group_by(cond, type) %>% 
  summarize(Counts=sum(N),
            Condition=unique(Condition),
            Celltype=unique(Celltype)) %>%
  as.data.frame()

#combine the three mitoribo_agg tables
mitoribo_agg <- bind_rows(mitoribo_agg_HEK_ex,mitoribo_agg_PBMC_ex,mitoribo_agg_Tissue_ex)

mitoribo_agg_inex <- bind_rows(mitoribo_agg_HEK_inex,mitoribo_agg_PBMC_inex,mitoribo_agg_Tissue_inex)


#add mito and ribo amounts to df
mapped <- mapped %>% 
  bind_rows(mitoribo_agg_inex) %>% 
  mutate(Condition=factor(Condition,levels=c( "Magnetic Beads",
                                               "Incubation",
                                               "No Incubation",
                                               "Column")),
         type_map=if_else(as.character(type)%in% c("Coding","lncRNA","Mitochondrial","Ribosomal","Intron","User"),
                          "Intragenic",
                          type),
         type=factor(type,
                     levels=rev(c("Coding","Intron","lncRNA","Mitochondrial","Ribosomal","Intergenic","Ambiguity","User")))) %>% 
  filter(type != "User")


mapped_type<-mapped %>% 
  group_by(Condition,Celltype,type_map,cond) %>% 
  summarize(Counts=sum(Counts))
```

#### 9.2 Make Plots

``` r

feat_cols<-c("#F08C4B", "#E5DFCC", "#3A8DB5", "#556f44", "#9dc183", "#4D8C57", "#256EA0","dodgerblue4")

names(feat_cols)<-c("Ambiguity","Intergenic","Intron","Ribosomal","Mitochondrial","lncRNA","Exon","Intragenic")

#Figure 2a
plot_type_map <- ggplot(subset(mapped_type, Condition %in% c("Magnetic Beads", "Column")), aes(x=Condition, y=Counts, fill=type_map))+
    geom_bar(stat="identity", position = "fill")+
    facet_grid(Celltype~., scales="free") +
    ylab("Fraction of Assigned Reads") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    coord_flip() +
    scale_fill_manual(values = feat_cols)+
    guides(fill = guide_legend(reverse=T))+
    theme(legend.title = element_blank(), 
          legend.position="bottom", 
          axis.title.y = element_blank(), 
          #axis.text.y = element_blank(),
          strip.text.y = element_text(angle = 360),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank()) 
plot_type_map
```

![](Lysis_features_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r

#Figure 2b
mapped_gtype<-mapped %>% 
  filter(type %in% c("Coding","lncRNA","Mitochondrial","Ribosomal","Intron")) %>%
  mutate(type = as.character(type), type = if_else(type == "Coding", "Exon", type), type=factor(type,
  levels=rev(c("Exon","Intron","lncRNA","Mitochondrial","Ribosomal","Intergenic","Ambiguity"))))

plot_gtype_map <- ggplot(subset(mapped_gtype, Condition %in% c("Magnetic Beads", "Column")), aes(x=Condition, y=Counts, fill=type))+
    geom_bar(stat="identity", position = "fill")+
    facet_grid(Celltype~., scales="free") +
    ylab("Fraction of Assigned Reads") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    coord_flip() +
    scale_fill_manual(values = feat_cols)+
    guides(fill = guide_legend(reverse=T))+
    theme(legend.title = element_blank(), 
          legend.position="bottom", 
          axis.title.y = element_blank(), 
          #axis.text.y = element_blank(),
          strip.text.y = element_text(angle = 360),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank()) 
plot_gtype_map
```

![](Lysis_features_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r

mapped_HEK<-mapped %>% 
  mutate(type = as.character(type), type = if_else(type == "Coding", "Exon", type), type=factor(type,
  levels=rev(c("Exon","Intron","lncRNA","Mitochondrial","Ribosomal","Intergenic","Ambiguity"))))


#Figure Supp_feat_HEK_all
plot_feat_map_HEK <- ggplot(subset(mapped_HEK, Celltype == "HEK"), aes(x=Condition, y=Counts, fill=type))+
    geom_bar(stat="identity", position = "fill")+
    facet_grid(Celltype~., scales="free") +
    ylab("Fraction of Assigned Reads") +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    coord_flip() +
    scale_fill_manual(values = feat_cols)+
    guides(fill = guide_legend(reverse=T))+
    theme_pub+
    theme(legend.title = element_blank(), 
          legend.position="bottom", 
          axis.title.y = element_blank(), 
          #axis.text.y = element_blank(),
          strip.text.y = element_text(angle = 360),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank()) 

plot_feat_map_HEK
```

![](Lysis_features_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->
