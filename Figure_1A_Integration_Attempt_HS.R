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

> sessionInfo()
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] metap_1.9                   paletteer_1.5.0             ggsci_3.0.0                 pals_1.8                   
 [5] RColorBrewer_1.1-3          plyr_1.8.9                  reshape2_1.4.4              stringr_1.5.1              
 [9] hrbrthemes_0.8.0            survminer_0.4.9             survival_3.5-5              irlba_2.3.5.1              
[13] rliger_1.0.1                ggpubr_0.6.0                pheatmap_1.0.12             future_1.33.0              
[17] clustree_0.5.1              ggraph_2.1.0                textshape_1.7.3             tidyr_1.3.0                
[21] ggrepel_0.9.4               data.table_1.14.8           ggplot2_3.4.4               dplyr_1.1.4                
[25] cowplot_1.1.1               Seurat_5.0.1                SeuratObject_5.0.1          sp_2.1-2                   
[29] batchelor_1.16.0            SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2 Biobase_2.60.0             
[33] GenomicRanges_1.52.1        GenomeInfoDb_1.36.4         IRanges_2.34.1              S4Vectors_0.38.2           
[37] BiocGenerics_0.46.0         MatrixGenerics_1.12.3       matrixStats_1.1.0           Matrix_1.6-4               
[41] SeuratWrappers_0.3.2       

loaded via a namespace (and not attached):
  [1] spatstat.sparse_3.0-3     bitops_1.0-7              httr_1.4.7                doParallel_1.0.17        
  [5] numDeriv_2016.8-1.1       tools_4.3.1               sctransform_0.4.1         backports_1.4.1          
  [9] utf8_1.2.4                R6_2.5.1                  ResidualMatrix_1.10.0     sn_2.1.1                 
 [13] lazyeval_0.2.2            uwot_0.1.16               withr_2.5.2               prettyunits_1.2.0        
 [17] gridExtra_2.3             progressr_0.14.0          textshaping_0.3.7         cli_3.6.1                
 [21] spatstat.explore_3.2-5    fastDummies_1.7.3         sandwich_3.0-2            labeling_0.4.3           
 [25] mvtnorm_1.2-4             survMisc_0.5.6            spatstat.data_3.0-3       ggridges_0.5.4           
 [29] pbapply_1.7-2             systemfonts_1.0.5         gfonts_0.2.0              R.utils_2.12.3           
 [33] dichromat_2.0-0.1         parallelly_1.36.0         plotrix_3.8-4             limma_3.56.2             
 [37] maps_3.4.1.1              rstudioapi_0.15.0         httpcode_0.3.0            FNN_1.1.3.2              
 [41] generics_0.1.3            ica_1.0-3                 spatstat.random_3.2-2     car_3.1-2                
 [45] fansi_1.0.5               abind_1.4-5               R.methodsS3_1.8.2         lifecycle_1.0.4          
 [49] multcomp_1.4-25           yaml_2.3.7                carData_3.0-5             mathjaxr_1.6-0           
 [53] Rtsne_0.16                grid_4.3.1                promises_1.2.1            crayon_1.5.2             
 [57] miniUI_0.1.1.1            lattice_0.21-8            beachmat_2.16.0           mapproj_1.2.11           
 [61] pillar_1.9.0              knitr_1.45                future.apply_1.11.0       codetools_0.2-19         
 [65] leiden_0.4.3.1            mutoss_0.1-13             glue_1.6.2                fontLiberation_0.1.0     
 [69] remotes_2.4.2.1           Rdpack_2.6                vctrs_0.6.5               png_0.1-8                
 [73] spam_2.10-0               gtable_0.3.4              rematch2_2.1.2            xfun_0.41                
 [77] rbibutils_2.2.16          S4Arrays_1.0.6            mime_0.12                 tidygraph_1.2.3          
 [81] iterators_1.0.14          KMsurv_0.1-5              TH.data_1.1-2             ellipsis_0.3.2           
 [85] fitdistrplus_1.1-11       ROCR_1.0-11               nlme_3.1-162              fontquiver_0.2.1         
 [89] bit64_4.0.5               RcppAnnoy_0.0.21          rprojroot_2.0.4           KernSmooth_2.23-21       
 [93] colorspace_2.1-0          mnormt_2.1.1              tidyselect_1.2.0          processx_3.8.2           
 [97] extrafontdb_1.0           bit_4.0.5                 compiler_4.3.1            curl_5.1.0               
[101] BiocNeighbors_1.18.0      hdf5r_1.3.8               TFisher_0.2.0             fontBitstreamVera_0.1.1  
[105] desc_1.4.2                DelayedArray_0.26.7       plotly_4.10.3             scales_1.3.0             
[109] lmtest_0.9-40             callr_3.7.3               digest_0.6.33             goftest_1.2-3            
[113] spatstat.utils_3.0-4      rmarkdown_2.25            XVector_0.40.0            htmltools_0.5.7          
[117] pkgconfig_2.0.3           extrafont_0.19            sparseMatrixStats_1.12.2  fastmap_1.1.1            
[121] rlang_1.1.2               htmlwidgets_1.6.4         shiny_1.8.0               DelayedMatrixStats_1.22.6
[125] farver_2.1.1              zoo_1.8-12                jsonlite_1.8.8            BiocParallel_1.34.2      
[129] mclust_6.0.1              R.oo_1.25.0               BiocSingular_1.16.0       RCurl_1.98-1.12          
[133] magrittr_2.0.3            scuttle_1.10.3            GenomeInfoDbData_1.2.10   dotCall64_1.1-1          
[137] patchwork_1.1.3           munsell_0.5.0             Rcpp_1.0.11               gdtools_0.3.4            
[141] viridis_0.6.4             reticulate_1.34.0         stringi_1.8.2             zlibbioc_1.46.0          
[145] MASS_7.3-60               pkgbuild_1.4.2            parallel_4.3.1            listenv_0.9.0            
[149] deldir_2.0-2              graphlayouts_1.0.2        splines_4.3.1             multtest_2.56.0          
[153] tensor_1.5                qqconf_1.3.2              ps_1.7.5                  igraph_1.5.1             
[157] spatstat.geom_3.2-7       ggsignif_0.6.4            RcppHNSW_0.5.0            ScaledMatrix_1.8.1       
[161] crul_1.4.0                evaluate_0.23             BiocManager_1.30.22       foreach_1.5.2            
[165] tweenr_2.0.2              httpuv_1.6.13             Rttf2pt1_1.3.12           RANN_2.6.1               
[169] purrr_1.0.2               polyclip_1.10-6           km.ci_0.5-6               scattermore_1.2          
[173] ggforce_0.4.1             rsvd_1.0.5                broom_1.0.5               xtable_1.8-4             
[177] RSpectra_0.16-1           rstatix_0.7.2             later_1.3.2               ragg_1.2.6               
[181] viridisLite_0.4.2         tibble_3.2.1              cluster_2.1.4             globals_0.16.2  
