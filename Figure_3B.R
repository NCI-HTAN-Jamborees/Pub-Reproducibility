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

###Load the SCLC Epithelial cell object
###This dataset can be found here: https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec

####Figure 3B Barchart of cluster by donor
dge <- readRDS("./Download_SCLC_Epithelial.rds")

library(reshape2)
dfr <- table(dge$clusters, dge$donor_id)
write.table(dfr, "SCLC_composition_bycluster.txt", sep="\t", col.names = T, row.names = T)
dfr <- read.table("./SCLC_composition_bycluster.txt", sep = "\t", header = T, row.names = 1)
dfr <- as.data.frame((dfr)) 
dfr$category <- row.names(dfr)
mdfr <- melt(dfr, id.vars = "category")
library(scales)
(p <- ggplot(mdfr, aes(category, value, fill = variable)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_y_continuous(labels = percent) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + RotatedAxis()
)
