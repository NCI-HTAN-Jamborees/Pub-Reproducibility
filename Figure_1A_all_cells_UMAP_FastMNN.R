#load the data which is downloaded from the following link
#https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec
#the "Combined samples" dataset in RDS format, and renamed it as "combined_data.rds"

#this code is adapted from the original file from the author's github
#https://github.com/dpeerlab/SCLC_atlas-HTAN/blob/main/scripts/run_fastmnn.R

##Zheng Xia, OHSU, 12/06/2023

rm(list=ls())
suppressMessages(library(Matrix))
suppressMessages(library(batchelor))
suppressMessages(library(BiocParallel))
suppressMessages(library(BiocParallel))
suppressMessages(library(Seurat))
ncores = 16
SnowParam(workers = ncores)
suppressMessages(library(irlba))
library(dplyr)


combined_data = readRDS("combined_data.rds")

num_hvg <- 5000 # 
knn <- 30

combined_data <- NormalizeData(combined_data)
combined_data <- FindVariableFeatures(combined_data, selection.method = "vst", nfeatures = num_hvg)
hvg = VariableFeatures(combined_data)



batch = combined_data$batch

combined_data <- NormalizeData(combined_data, normalization.method = "RC")
counts0=combined_data@assays$RNA@data


obs_df = combined_data@meta.data

obs_df = obs_df %>% dplyr::select(batch, HTAN_Participant_ID, histo, tissue)
obs_df2 = obs_df[!duplicated(obs_df$batch),] %>% dplyr::arrange(batch)

num_tissue = sapply(sort(unique(obs_df2$HTAN_Participant_ID)), function(x) length(unique(obs_df2[obs_df2$HTAN_Participant_ID==x, 'tissue'])))
num_cell = table(obs_df$batch)

obs_df$num_tissue = num_tissue[obs_df$HTAN_Participant_ID]
obs_df$num_cell = num_cell[obs_df$batch]

meta = obs_df %>% dplyr::rename(sample=batch)                


order_df = meta[!duplicated(meta$sample), ]
order_df$histo = gsub('normal','LUAD',order_df$histo) #RELABEL NORMAL TO LUAD SINCE THEY ARE NORMAL LUNG ADJACENT TO LUAD

#HIERARCHICAL MERGING STRATEGY: MERGE FIRST BY PATIENT THEN BY HISTOLOGY. MERGE PATIENTS WITH MULTIPLE TISSUE SITES FIRST. WITHIN PATIENT, MERGE SAMPLES WITH HIGHEST NUMBER OF CELLS FIRST
order_df = order_df[order(order_df$histo, order_df$num_tissue, order_df$HTAN_Participant_ID, order_df$num_cell, decreasing = TRUE),]

batch_order = sort(unique(order_df$sample))

merge.order = list()
for (i in unique(order_df$HTAN_Participant_ID)) {
  merge.order[[i]] = which(order_df$HTAN_Participant_ID==i)
}

merge.order2 = list()
for (i in unique(order_df$histo)) {
  tmp = order_df$HTAN_Participant_ID[order_df$histo==i]
  merge.order2[[i]] = merge.order[tmp[!duplicated(tmp)]]
}

counts = log2(counts0[hvg,] + 0.1)

batch = gsub('_[0-9]+$','',colnames(counts))

counts_list = lapply(split(seq_along(batch), batch), function(m, ind) m[,ind], m=counts)[order(unique(batch))]

counts_list = counts_list[batch_order]

#PLEASE NOTE CAN TOGGLE COSINE NORMALIZATION TO ON, BUT THIS WAS NOT USED FOR THE SCLC ATLAS 
set.seed(2023)
out = fastMNN(counts_list, batch = batch_order, k=knn, merge.order = merge.order, cos.norm=F, d = 50, get.variance=T)
save(out, file="fastmnn_out2.Rdata")


#load("fastmnn_out2.Rdata")
correct = reducedDims(out)$corrected
correct = correct[match(colnames(combined_data), rownames(correct)),]

colnames(correct) <- paste0('fastMNN_',1:50) 
combined_data[["fastMNN"]] <- CreateDimReducObject(embeddings=as.matrix(correct), assay='RNA',  key = "fastMNN_") 


combined_data = FindNeighbors(combined_data, reduction ='fastMNN', dims = 1:50)
combined_data <- FindClusters(combined_data, resolution = 0.5)

combined_data <- RunUMAP(object = combined_data, dims = 1:50, reduction = 'fastMNN')
plot_1 <- DimPlot(combined_data,reduction = 'umap' , label = T, group.by = "cell_type_coarse" ,repel = TRUE, 
                   label.box = TRUE, label.size=5)



#umap based on the pca from the rds object from the original paper
combined_data = readRDS("combined_data.rds")
combined_data <- RunUMAP(object = combined_data, dims = 1:50, reduction = 'pca')
DimPlot(combined_data,reduction = 'umap' , label = T, group.by = "cell_type_coarse" ,repel = TRUE, 
        label.box = TRUE, label.size=5)

##done

# > sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.3                 irlba_2.3.5.1               SeuratObject_4.1.4          Seurat_4.4.0               
# [5] BiocParallel_1.34.2         batchelor_1.16.0            SingleCellExperiment_1.22.0 SummarizedExperiment_1.30.2
# [9] Biobase_2.60.0              GenomicRanges_1.52.1        GenomeInfoDb_1.36.4         IRanges_2.34.1             
# [13] S4Vectors_0.38.2            BiocGenerics_0.46.0         MatrixGenerics_1.12.3       matrixStats_1.0.0          
# [17] Matrix_1.6-1.1             
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3        jsonlite_1.8.7            magrittr_2.0.3            spatstat.utils_3.0-3      farver_2.1.1             
# [6] rmarkdown_2.25            zlibbioc_1.46.0           vctrs_0.6.4               ROCR_1.0-11               spatstat.explore_3.2-3   
# [11] DelayedMatrixStats_1.22.6 RCurl_1.98-1.12           htmltools_0.5.6.1         S4Arrays_1.0.6            BiocNeighbors_1.18.0     
# [16] sctransform_0.4.0         parallelly_1.36.0         KernSmooth_2.23-21        htmlwidgets_1.6.2         ica_1.0-3                
# [21] plyr_1.8.9                plotly_4.10.2             zoo_1.8-12                ResidualMatrix_1.10.0     igraph_1.5.1             
# [26] mime_0.12                 lifecycle_1.0.3           pkgconfig_2.0.3           rsvd_1.0.5                R6_2.5.1                 
# [31] fastmap_1.1.1             GenomeInfoDbData_1.2.10   fitdistrplus_1.1-11       future_1.33.0             shiny_1.7.5.1            
# [36] digest_0.6.33             colorspace_2.1-0          patchwork_1.1.3           tensor_1.5                beachmat_2.16.0          
# [41] labeling_0.4.3            progressr_0.14.0          spatstat.sparse_3.0-2     fansi_1.0.5               polyclip_1.10-6          
# [46] httr_1.4.7                abind_1.4-5               compiler_4.3.1            withr_2.5.1               MASS_7.3-60              
# [51] DelayedArray_0.26.7       tools_4.3.1               lmtest_0.9-40             httpuv_1.6.11             future.apply_1.11.0      
# [56] goftest_1.2-3             glue_1.6.2                nlme_3.1-162              promises_1.2.1            grid_4.3.1               
# [61] Rtsne_0.16                cluster_2.1.4             reshape2_1.4.4            generics_0.1.3            spatstat.data_3.0-1      
# [66] gtable_0.3.4              tidyr_1.3.0               data.table_1.14.8         BiocSingular_1.16.0       ScaledMatrix_1.8.1       
# [71] sp_2.1-1                  utf8_1.2.3                XVector_0.40.0            spatstat.geom_3.2-5       RcppAnnoy_0.0.21         
# [76] ggrepel_0.9.4             RANN_2.6.1                pillar_1.9.0              stringr_1.5.0             later_1.3.1              
# [81] splines_4.3.1             lattice_0.21-8            deldir_1.0-9              survival_3.5-5            tidyselect_1.2.0         
# [86] miniUI_0.1.1.1            scuttle_1.10.3            pbapply_1.7-2             knitr_1.44                gridExtra_2.3            
# [91] scattermore_1.2           xfun_0.40                 stringi_1.7.12            lazyeval_0.2.2            yaml_2.3.7               
# [96] evaluate_0.22             codetools_0.2-19          tibble_3.2.1              cli_3.6.1                 uwot_0.1.16              
# [101] xtable_1.8-4              reticulate_1.34.0         munsell_0.5.0             Rcpp_1.0.11               spatstat.random_3.1-6    
# [106] globals_0.16.2            png_0.1-8                 parallel_4.3.1            ellipsis_0.3.2            ggplot2_3.4.4            
# [111] sparseMatrixStats_1.12.2  bitops_1.0-7              listenv_0.9.0             viridisLite_0.4.2         scales_1.2.1             
# [116] ggridges_0.5.4            leiden_0.4.3              purrr_1.0.2               crayon_1.5.2              rlang_1.1.1              
# [121] cowplot_1.1.1 

