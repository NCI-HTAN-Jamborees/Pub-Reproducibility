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

####Figure S2A
###Subtype Diversity and Frequency
dge<-readRDS("./SCLC_Epithelial.rds")

####Figure 1H
p1<-FeaturePlot(dge,features=c("ASCL1"),cols = c("purple","red","yellow"))
p2<-FeaturePlot(dge,features=c("NEUROD1"),cols = c("purple","red","yellow"))
p3<-FeaturePlot(dge,features=c("POU2F3"),cols = c("purple","red","yellow"))
p4<-FeaturePlot(dge,features=c("YAP1"),cols = c("purple","red","yellow"))
grid.arrange(p1, p2, p3, p4, col = 2)
####Done
