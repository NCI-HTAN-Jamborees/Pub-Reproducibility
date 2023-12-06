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

# > sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: x86_64-apple-darwin20 (64-bit)
# Running under: macOS Big Sur ... 10.16
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scales_1.2.1       mgcv_1.8-42        nlme_3.1-162       ggtern_3.4.2       pals_1.7          
# [6] RColorBrewer_1.1-3 entropy_1.3.1      plyr_1.8.8         reshape2_1.4.4     stringr_1.5.0     
# [11] jjb_0.1.1          ggpubr_0.6.0       pheatmap_1.0.12    future_1.33.0      clustree_0.5.0    
# [16] ggraph_2.1.0       textshape_1.7.3    tidyr_1.3.0        ggrepel_0.9.3      data.table_1.14.8 
# [21] ggplot2_3.4.2      dplyr_1.1.2        cowplot_1.1.1      SeuratObject_4.1.3 Seurat_4.3.0.1    
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                    matrixStats_1.0.0           spatstat.sparse_3.0-2      
# [4] bitops_1.0-7                httr_1.4.6                  latex2exp_0.9.6            
# [7] tools_4.3.1                 sctransform_0.3.5           backports_1.4.1            
# [10] utf8_1.2.3                  R6_2.5.1                    lazyeval_0.2.2             
# [13] uwot_0.1.16                 withr_2.5.0                 sp_2.0-0                   
# [16] gridExtra_2.3               bayesm_3.1-6                progressr_0.13.0           
# [19] cli_3.6.1                   Biobase_2.60.0              textshaping_0.3.6          
# [22] spatstat.explore_3.2-1      officer_0.6.2               labeling_0.4.2             
# [25] robustbase_0.99-0           spatstat.data_3.0-1         ggridges_0.5.4             
# [28] pbapply_1.7-2               askpass_1.1                 systemfonts_1.0.4          
# [31] gfonts_0.2.0                dichromat_2.0-0.1           scater_1.28.0              
# [34] parallelly_1.36.0           maps_3.4.1                  rstudioapi_0.15.0          
# [37] httpcode_0.3.0              generics_0.1.3              ica_1.0-3                  
# [40] spatstat.random_3.1-5       car_3.1-2                   zip_2.3.0                  
# [43] Matrix_1.6-0                ggbeeswarm_0.7.2            fansi_1.0.4                
# [46] S4Vectors_0.38.1            abind_1.4-5                 lifecycle_1.0.3            
# [49] yaml_2.3.7                  carData_3.0-5               SummarizedExperiment_1.30.2
# [52] Rtsne_0.16                  promises_1.2.0.1            crayon_1.5.2               
# [55] miniUI_0.1.1.1              lattice_0.21-8              beachmat_2.16.0            
# [58] mapproj_1.2.11              pillar_1.9.0                knitr_1.43                 
# [61] GenomicRanges_1.52.0        future.apply_1.11.0         codetools_0.2-19           
# [64] leiden_0.4.3                glue_1.6.2                  fontLiberation_0.1.0       
# [67] vctrs_0.6.3                 png_0.1-8                   gtable_0.3.3               
# [70] xfun_0.39                   S4Arrays_1.0.5              mime_0.12                  
# [73] tidygraph_1.2.3             survival_3.5-5              SingleCellExperiment_1.22.0
# [76] ellipsis_0.3.2              fitdistrplus_1.1-11         ROCR_1.0-11                
# [79] usethis_2.2.2               fontquiver_0.2.1            RcppAnnoy_0.0.21           
# [82] GenomeInfoDb_1.36.1         tensorA_0.36.2              irlba_2.3.5.1              
# [85] vipor_0.4.5                 KernSmooth_2.23-21          colorspace_2.1-0           
# [88] BiocGenerics_0.46.0         tidyselect_1.2.0            compiler_4.3.1             
# [91] curl_5.0.1                  compositions_2.0-6          BiocNeighbors_1.18.0       
# [94] flextable_0.9.2             xml2_1.3.5                  fontBitstreamVera_0.1.1    
# [97] DelayedArray_0.26.7         plotly_4.10.2               DEoptimR_1.1-3             
# [100] lmtest_0.9-40               hexbin_1.28.3               digest_0.6.33              
# [103] goftest_1.2-3               spatstat.utils_3.0-3        rmarkdown_2.23             
# [106] XVector_0.40.0              htmltools_0.5.5             pkgconfig_2.0.3            
# [109] sparseMatrixStats_1.12.2    MatrixGenerics_1.12.3       fastmap_1.1.1              
# [112] rlang_1.1.1                 htmlwidgets_1.6.2           shiny_1.7.4.1              
# [115] DelayedMatrixStats_1.22.1   farver_2.1.1                zoo_1.8-12                 
# [118] jsonlite_1.8.7              BiocParallel_1.34.2         BiocSingular_1.16.0        
# [121] RCurl_1.98-1.12             magrittr_2.0.3              scuttle_1.10.1             
# [124] GenomeInfoDbData_1.2.10     patchwork_1.1.2             munsell_0.5.0              
# [127] Rcpp_1.0.11                 viridis_0.6.4               proto_1.0.0                
# [130] gdtools_0.3.3               reticulate_1.31             stringi_1.7.12             
# [133] zlibbioc_1.46.0             MASS_7.3-60                 parallel_4.3.1             
# [136] listenv_0.9.0               deldir_1.0-9                graphlayouts_1.0.0         
# [139] splines_4.3.1               tensor_1.5                  igraph_1.5.0.1             
# [142] uuid_1.1-0                  spatstat.geom_3.2-4         ggsignif_0.6.4             
# [145] stats4_4.3.1                ScaledMatrix_1.8.1          crul_1.4.0                 
# [148] evaluate_0.21               tweenr_2.0.2                httpuv_1.6.11              
# [151] RANN_2.6.1                  openssl_2.1.0               purrr_1.0.1                
# [154] polyclip_1.10-4             scattermore_1.2             ggforce_0.4.1              
# [157] rsvd_1.0.5                  broom_1.0.5                 xtable_1.8-4               
# [160] rstatix_0.7.2               later_1.3.1                 viridisLite_0.4.2          
# [163] ragg_1.2.5                  tibble_3.2.1                beeswarm_0.4.0             
# [166] IRanges_2.34.1              cluster_2.1.4               globals_0.16.2             
