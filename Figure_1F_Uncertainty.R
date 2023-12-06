
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



#### Make Epithelial subset
### Load the all cell object
### This dataset can be found here: https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec

dge <- readRDS("./All_Cell.rds")
dge <- subset(dge,cells=colnames(dge)[dge$cell_type_coarse=="Epithelial"])

### This csv file is the output for the python code.
Train_data <- read.table("./Epithelial_metadata.csv",sep=",",header = T,row.names = 1)
dge@meta.data <- cbind(dge@meta.data,Train_data[,43:46])


####Uncertainty analysis Figure 1F
metadata_df = dge@meta.data
metadata_df$batch = as.character(metadata_df$batch)
metadata_df$batch[grepl('RU1108',metadata_df$batch)] = 'RU1108a'
metadata_df$subtype_uncertainty = apply(metadata_df[,c("pval_SCLC.A","pval_SCLC.N","pval_SCLC.P")], 1, function(x) entropy(x))
any(is.nan(metadata_df$subtype_uncertainty))
library(ggtern)

pval_df = metadata_df[,grepl('pval_SCLC.[ANP]',colnames(metadata_df))]
sigma = 0.02
pval_df = pval_df + sigma   #Allows jitter for better visualization
colnames(pval_df) = gsub('pval_','',colnames(pval_df))

pval_df = pval_df + cbind(rnorm(n = nrow(pval_df), mean = 0, sd = sigma), 
                          rnorm(n = nrow(pval_df), mean = 0, sd = sigma),
                          rnorm(n = nrow(pval_df), mean = 0, sd = sigma)) #Create jitter 

pval_df = pval_df - apply(pval_df, 1, function(x) min(c(x,0)) )
pval_df = pval_df/rowSums(pval_df)

col = rgb(pval_df$SCLC.A, pval_df$SCLC.N, pval_df$SCLC.P)

pval_gg = as.data.frame(cbind(pval_df, col))

gplot <- ggtern(data = pval_gg, 
                aes(y = SCLC.A, x = SCLC.N, z = SCLC.P, col=col),
                aes(x,y,z)) +
  geom_point(aes_string(col = 'col'), size=0.01) +
  labs(y='SCLC-A', x='SCLC-N', z='SCLC-P') + 
  scale_color_identity() + 
  theme_bw() +
  theme(tern.panel.mask.show = FALSE,
        axis.title = element_text(size=21),
        tern.axis.title.L = element_text(hjust = 0,vjust=1.8),
        tern.axis.title.R = element_text(hjust = 1,vjust=1.8),
        text = element_text(size=17))
gplot

####Done

