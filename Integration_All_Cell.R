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

#### Read in the merged processed object
#### This object is the result of Parse_Overall_Merge.R
#### For integration, use 50 PCs (same as the Chan et al. paper)
seurat<-readRDS("./Overall_HTAN.rds")

####LIGER integration by batch
####
library(rliger)
library(SeuratWrappers)

seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat, split.by = "batch", do.center = FALSE)
seurat <- RunOptimizeALS(seurat, k = 50, lambda = 10, split.by = "batch")
seurat <- RunQuantileNorm(seurat, split.by = "batch")
# You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
# `RunQuantileNorm` according to your needs
seurat <- FindNeighbors(seurat, reduction = "iNMF", dims = 1:ncol(Embeddings(seurat, "iNMF")))
seurat <- FindClusters(seurat, resolution = 0.3)
# Dimensional reduction and plotting
seurat <- RunUMAP(seurat, dims = 1:ncol(seurat[["iNMF"]]), reduction = "iNMF")
saveRDS(seurat,"./Integrated_LIGER.rds")

####Test the visualization for coarse cell type and fine cell type
DimPlot(seurat, reduction = "umap",label=T,group.by = "cell_type_coarse",raster = T)+NoLegend()
ggsave(file="Umap_raw_type_LIGER.pdf",width = 20,height = 20,units = "cm")
DimPlot(seurat, reduction = "umap",group.by="cell_type_fine",label=T,raster = T)
ggsave(file="Umap_fineID_LIGER.pdf",width = 20,height = 15)
DimPlot(seurat, reduction = "umap",group.by="donor_id",label=T,raster = T)
ggsave(file="Umap_donor_LIGER.pdf",width = 20,height = 15)

saveRDS(seurat,"./Integrated_LIGER.rds")
#### LIGER done

####Second, test MNN integration by batch
####Re-load the merged object
seurat<-readRDS("./Overall_HTAN.rds")
seurat_samples <- SplitObject(seurat, "batch")
seurat_mnn <- RunFastMNN(seurat_samples)
seurat[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn, "mnn")[colnames(seurat),], key="MNN_")
seurat <- RunUMAP(seurat, reduction = "mnn",dims = 1:50)

####Test the visualization for coarse cell type and fine cell type
DimPlot(seurat, reduction = "umap",label=T,group.by = "cell_type_coarse",raster = T)+NoLegend()
ggsave(file="Umap_raw_type_MNN.pdf",width = 20,height = 20,units = "cm")
DimPlot(seurat, reduction = "umap",group.by="cell_type_fine",label=T,raster = T)
ggsave(file="Umap_fineID_MNN.pdf",width = 20,height = 15)
DimPlot(seurat, reduction = "umap",group.by="donor_id",label=T,raster = T)
ggsave(file="Umap_donor_MNN.pdf",width = 20,height = 15)

saveRDS(seurat,"./Integrated_MNN.rds")
#### MNN done

