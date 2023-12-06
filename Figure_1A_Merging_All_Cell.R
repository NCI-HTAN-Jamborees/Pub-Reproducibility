
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

### Load the all cell object
### This dataset can be found here: https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec

#### Convert the geneID (rownames(dge)) to gene symbols using the meta features
#### Re-make the RNA assay
RNA <- dge@assays$RNA
RNA@counts@Dimnames[[1]] <- as.vector(Features$feature_name)
RNA@data@Dimnames[[1]] <- as.vector(Features$feature_name)
dge@assays$RNA <- RNA

#### Re-create a Seurat object for re-processing. 
#### Note: If users want to skip this and re-use the computed UMAP in the manuscript, JUMP to LINE 93

#### Start Re-processing 
dge<-CreateSeuratObject(counts = dge@assays$RNA@counts,meta.data = dge@meta.data)
#dge<-readRDS("/wynton/group/fhuang/Atlas/RPCA_Atlas_integrated_clustered.rds")
####Top_tier annotation
options(future.globals.maxSize = 80000 * 1024^2)
setwd("/Volumes/Backup_Plus/HTAN/")


#####QC
VlnPlot(object = dge, features = c("ngenes", "mito_frac", "RBP_frac"),group.by = "donor_id",ncol=1,pt.size = 0)
ggsave(file="QC_bydonor.pdf",width = 15,height = 10)

VlnPlot(object = dge, features = c("ngenes", "mito_frac", "RBP_frac"),group.by = "procedure",ncol=1,pt.size = 0)
ggsave(file="QC_byprocedure.pdf",width = 5,height = 8)

dge[["percent.mt"]] <- PercentageFeatureSet(object = dge, pattern = "^MT-")
dge[["percent.rps"]] <- PercentageFeatureSet(object = dge, pattern = "^RPS")
dge[["percent.rpl"]] <- PercentageFeatureSet(object = dge, pattern = "^RPL")
dge[["percent.ribosomal"]] <- dge[["percent.rps"]]+dge[["percent.rpl"]]


VlnPlot(object = dge, features = c("nFeature_RNA", "nCount_RNA"),ncol = 1,group.by = "donor_id",pt.size = 0)
ggsave(file="QC_bydonor_raw.pdf",width = 15,height = 8)

dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)

dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures =3000)

# all.genes <- rownames(x = dge)
dge <- ScaleData(object=dge,features=VariableFeatures(object = dge))
dge <- RunPCA(object = dge, features = VariableFeatures(object=dge),npcs = 100)


# Create a vector of gene names for downstream analysis
#gene.names.vector <- as.character("gene names")
# Run PCA on selected genes
#dge <- RunPCA(object = dge, features = gene.names.vector, do.print = TRUE, pcs.print = 1:5, 
#                 genes.print = 5)

slot(dge[["pca"]], "misc")
print(x = dge[["pca"]], dims = 1:5, nfeatures = 5)
mat <- Seurat::GetAssayData(dge,assay="RNA",slot="scale.data")
pca <-dge[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)

DimPlot(object = dge, reduction = "pca",group.by = "donor_id",raster = T)
ggsave(file="PCA_bydonor.pdf",width = 20,height = 20,units = "cm")
png("Elbowplot_dge.png")
ElbowPlot(dge,ndims = 100)
dev.off()
#### Use 100 PCs for UMAP (from the elbowplot)
n_pc=100
dge <- RunUMAP(dge, dims = 1:100)
#### Done re-processing and re-compute UMAP
#### Skip to here if want to use the paper-computed UMAP
DimPlot(dge, reduction = "umap",label=T,group.by = "cell_type_coarse",raster = T)+NoLegend()
ggsave(file="Umap_raw_type.pdf",width = 20,height = 20,units = "cm")


DimPlot(dge, reduction = "umap",group.by="cell_type_fine",label=T,raster = T)
ggsave(file="Umap_fineID.pdf",width = 20,height = 15)

DimPlot(dge, reduction = "umap",group.by="donor_id",label=T,raster = T)
ggsave(file="Umap_donor.pdf",width = 20,height = 15)


