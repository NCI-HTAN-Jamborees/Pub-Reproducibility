library(Seurat)
library(stringr)
library("ggplot2")
library("dplyr")
library('reshape2')
library(tidyr)
library(plyr)
library(ggrepel)
library(RColorBrewer)
library(pals)
library("ggsci")
library(paletteer)
library(metap)

seurat_data <- readRDS("./SCLC_epithelial_cells.rds")

#### switch gene annotation ####
df_gene_anno <- seurat_data@assays$RNA@meta.features
unique(rownames(seurat_data@assays$RNA@counts) == row.names(df_gene_anno))
rownames(seurat_data@assays$RNA@counts) = df_gene_anno$feature_name
rownames(seurat_data@assays$RNA@data) = df_gene_anno$feature_name

##### samples from paper
sample_interest = c(
  'RU1065C',
  'RU1195A',
  'RU779D',
  'RU1124A_LN',
  'RU1181B',
  'RU426B',
  'RU1215' )

#### perform differential expression between the recurrent cluster and the others
# Idents(object = seurat_data) <- "batch"
# DimPlot(seurat_data, reduction = "umap",label = T)
# DimPlot(subset(x = seurat_data, subset = batch == "RU1215"), reduction = "umap",label = T)

Idents(object = seurat_data) <- "clusters"
DimPlot(seurat_data, reduction = "umap",label = T)


df_deg <- FindMarkers(object = subset(x = seurat_data, subset = batch == sample_interest[1]),
                      slot = "data",
                      ident.1 = 22)
df_deg_use <- data.frame(Gene=row.names(df_deg),Value=df_deg$p_val_adj)
colnames(df_deg_use) <- c("Gene",sample_interest[1])
df_deg_com <-df_deg_use
for ( sample in sample_interest[-1] ){
  print(paste0("DEG analysis for ",sample))
  
  df_deg <- FindMarkers(object = subset(x = seurat_data, subset = batch == sample),
                        slot = "data",
                        ident.1 = 22)
  df_deg_use <- data.frame(Gene=row.names(df_deg),Value=df_deg$p_val_adj)
  colnames(df_deg_use) <- c("Gene",sample)
  df_deg_com <- left_join(df_deg_com,df_deg_use)
}

df_all = na.omit(df_deg_com)  
#df_deg_com[is.na(df_deg_com)] = 1
#df_all=df_deg_com
rownames(df_all) = df_all$Gene
df_all = df_all %>% dplyr::select(-Gene)

### The adjusted p values for DE within each tumor are combined using Edgington’s method
comb_pval = apply(df_all, 1, function(x) sump(x)$p)
comb_pval = sort(comb_pval)

plot_df = data.frame(gene = names(comb_pval), pval = comb_pval)

plot_df$rank = 1:nrow(plot_df) 
plot_df$log_pval = -log2(plot_df$pval)
plot_df$top = plot_df$rank != 1
plot_df$top_label = plot_df$gene
plot_df$top_label[20:nrow(plot_df)] = ''

p = plot_df %>% ggplot(aes(x = rank,y = log_pval, label = top_label)) + geom_point(aes(color = top))  +
  theme_bw() + xlab('Rank') + ylab('log2 (combined p-value)')  +
  geom_text_repel(segment.color = 'lightgray', segment.alpha= 0.5, max.overlaps = 100) +
  theme(legend.position = 'none',panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggsave("Fig3F_Genes_ordered_by_recurrence_score.pdf",p,height = 4.3,width = 5.3)


# > sessionInfo()
# R version 4.3.2 (2023-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.6.2
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Los_Angeles
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] metap_1.9           BiocManager_1.30.22 paletteer_1.5.0     ggsci_3.0.0         pals_1.8            RColorBrewer_1.1-3  ggrepel_0.9.4      
# [8] plyr_1.8.9          tidyr_1.3.0         reshape2_1.4.4      dplyr_1.1.4         stringr_1.5.1       Seurat_5.0.1        SeuratObject_5.0.1 
# [15] sp_2.1-2            hrbrthemes_0.8.0    survminer_0.4.9     ggpubr_0.6.0        ggplot2_3.4.4       survival_3.5-7     
# 
# loaded via a namespace (and not attached):
#   [1] ggtext_0.1.2            fs_1.6.3                matrixStats_1.1.0       spatstat.sparse_3.0-3   devtools_2.4.5          httr_1.4.7             
# [7] numDeriv_2016.8-1.1     profvis_0.3.8           tools_4.3.2             sctransform_0.4.1       backports_1.4.1         utf8_1.2.4             
# [13] R6_2.5.1                lazyeval_0.2.2          uwot_0.1.16             mgcv_1.9-0              sn_2.1.1                urlchecker_1.0.1       
# [19] withr_2.5.2             prettyunits_1.2.0       gridExtra_2.3           progressr_0.14.0        Biobase_2.62.0          cli_3.6.1              
# [25] textshaping_0.3.7       spatstat.explore_3.2-5  fastDummies_1.7.3       sandwich_3.0-2          labeling_0.4.3          prismatic_1.1.1        
# [31] mvtnorm_1.2-4           survMisc_0.5.6          spatstat.data_3.0-3     ggridges_0.5.4          pbapply_1.7-2           systemfonts_1.0.5      
# [37] commonmark_1.9.0        gfonts_0.2.0            dichromat_2.0-0.1       parallelly_1.36.0       sessioninfo_1.2.2       plotrix_3.8-4          
# [43] maps_3.4.1.1            rstudioapi_0.15.0       httpcode_0.3.0          generics_0.1.3          ica_1.0-3               spatstat.random_3.2-2  
# [49] car_3.1-2               Matrix_1.6-4            fansi_1.0.5             abind_1.4-5             lifecycle_1.0.4         multcomp_1.4-25        
# [55] yaml_2.3.7              carData_3.0-5           mathjaxr_1.6-0          Rtsne_0.16              grid_4.3.2              promises_1.2.1         
# [61] crayon_1.5.2            miniUI_0.1.1.1          lattice_0.22-5          cowplot_1.1.1           mapproj_1.2.11          pillar_1.9.0           
# [67] knitr_1.45              future.apply_1.11.0     codetools_0.2-19        leiden_0.4.3.1          mutoss_0.1-13           glue_1.6.2             
# [73] fontLiberation_0.1.0    data.table_1.14.8       remotes_2.4.2.1         vctrs_0.6.5             png_0.1-8               spam_2.10-0            
# [79] Rdpack_2.6              gtable_0.3.4            rematch2_2.1.2          cachem_1.0.8            xfun_0.41               rbibutils_2.2.16       
# [85] mime_0.12               KMsurv_0.1-5            ellipsis_0.3.2          fitdistrplus_1.1-11     TH.data_1.1-2           ROCR_1.0-11            
# [91] nlme_3.1-164            usethis_2.2.2           fontquiver_0.2.1        RcppAnnoy_0.0.21        irlba_2.3.5.1           KernSmooth_2.23-22     
# [97] BiocGenerics_0.48.1     colorspace_2.1-0        mnormt_2.1.1            tidyselect_1.2.0        processx_3.8.2          compiler_4.3.2         
# [103] extrafontdb_1.0         curl_5.1.0              xml2_1.3.5              TFisher_0.2.0           fontBitstreamVera_0.1.1 plotly_4.10.3          
# [109] scales_1.3.0            lmtest_0.9-40           callr_3.7.3             digest_0.6.33           goftest_1.2-3           presto_1.0.0           
# [115] spatstat.utils_3.0-4    rmarkdown_2.25          htmltools_0.5.7         pkgconfig_2.0.3         extrafont_0.19          fastmap_1.1.1          
# [121] rlang_1.1.2             htmlwidgets_1.6.3       shiny_1.8.0             farver_2.1.1            zoo_1.8-12              jsonlite_1.8.7         
# [127] magrittr_2.0.3          dotCall64_1.1-1         patchwork_1.1.3         munsell_0.5.0           Rcpp_1.0.11             gdtools_0.3.4          
# [133] reticulate_1.34.0       stringi_1.8.2           MASS_7.3-60             pkgbuild_1.4.2          parallel_4.3.2          listenv_0.9.0          
# [139] deldir_2.0-2            splines_4.3.2           multtest_2.58.0         gridtext_0.1.5          tensor_1.5              qqconf_1.3.2           
# [145] ps_1.7.5                igraph_1.5.1            spatstat.geom_3.2-7     ggsignif_0.6.4          markdown_1.11           RcppHNSW_0.5.0         
# [151] stats4_4.3.2            pkgload_1.3.3           crul_1.4.0              evaluate_0.23           httpuv_1.6.12           Rttf2pt1_1.3.12        
# [157] RANN_2.6.1              purrr_1.0.2             polyclip_1.10-6         future_1.33.0           km.ci_0.5-6             scattermore_1.2        
# [163] broom_1.0.5             xtable_1.8-4            RSpectra_0.16-1         rstatix_0.7.2           later_1.3.1             viridisLite_0.4.2      
# [169] ragg_1.2.6              tibble_3.2.1            memoise_2.0.1           cluster_2.1.6           globals_0.16.2 









