library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(tidyr)
library(textshape)
library(clustree)
library(future)
library(pheatmap)
library(ggpubr)
library(ggtern)

###Load the SCLC Epithelial cell object
###This dataset can be found here: https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec

####Figure 2A dotplot

dge<-readRDS("./SCLC_Epithelial.rds")

dge<-SetIdent(dge,value = as.vector(dge$SCLC_subtype))
dge@active.ident<-factor(dge@active.ident,levels = c("SCLC-A","SCLC-N","SCLC-P"))
Genelist<-c("ASCL1","TFF3","SOX4","KLF6","NEUROD1","NEUROD4","POU2F3","ASCL2",
            "SOX9","YBX3","HOXC10","YAP1","MYC","MYCN","MYCL","NOTCH1","HES6",
            "HES1","TCF4","MARCKS","ADCYAP1","NRXN1","SEMA6A","EFNB1","EPHB2",
            "NRP2","OLFM2","SSTR2","SST","SOX7","INHBA","HIF1A","VEGFA","FOXO3","VIM",
            "COL1A2","TWIST1","ZEB1","MMP2","MMP14","TGFB1","TGFBR1","TGFBR3","BMP3","BMP7",
            "BMPR1A","BMPR2","STAT3","IL6R","IL11RA","IL13RA1","TNF","PHLDA1","SMAD3",
            "SPHK1","NFAT5")

DotPlot(dge, features = rev(Genelist), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))+coord_flip()
ggsave(file="Figure2A_Replicate.pdf",width = 8,height = 15)
