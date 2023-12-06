
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

####Figure 2B scatter plot with spline
dge <- readRDS("./Download_SCLC_Epithelial.rds")

Features <- dge@assays$RNA@meta.features
RNA <- dge@assays$RNA
RNA@counts@Dimnames[[1]] <- as.vector(Features$feature_name)
RNA@data@Dimnames[[1]] <- as.vector(Features$feature_name)
dge@assays$RNA <- RNA

###Compute signature scores for these four gene sets
###These gene sets can be downloaded from msigdb

Featurename<-c("AXONOGENESIS","NEURON_DIFFERENTIATION","NEUROPEPTIDE_RECEPTOR","EMT")
Features<-list()
Features[[1]] <- read.table("./AXONOGENESIS.v2023.2.Hs.gmt",sep="\t",header = F)
Features[[1]] <- as.character(Features[[1]][1,])
Features[[2]] <- read.table("./GOBP_NEURON_DIFFERENTIATION.v2023.2.Hs.gmt",sep="\t",header = F)
Features[[2]] <- as.character(Features[[2]][1,])
Features[[3]] <- read.table("./GOMF_NEUROPEPTIDE_RECEPTOR_ACTIVITY.v2023.2.Hs.gmt",sep="\t",header = F)
Features[[3]] <- as.character(Features[[3]][1,])
Features[[4]] <- read.table("./HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2023.2.Hs.gmt",sep="\t",header = F)
Features[[4]] <- as.character(Features[[4]][1,])

for (i in 1:length(Features)) {
  gene.set <- Features[[i]][1:length(Features[[i]])]
  dge <- AddModuleScore(object = dge, features = list(gene.set), ctrl = 5, name = Featurename[i])
}

dge$ASCL1_Expression <- FetchData(dge, vars = "ASCL1")
dge$NEUROD_genes <- dge$`NEUROD genes`

scale_zscore <- function(DF) {
  DF_tmp <- DF - min(DF)
  DF_tmp <- DF_tmp / max(DF)
}

dge$ASCL1_Expression <- scale_zscore(dge$ASCL1_Expression)
dge$NEUROD_genes <- scale_zscore(dge$NEUROD_genes)
dge$AXONOGENESIS1 <- scale_zscore(dge$AXONOGENESIS1)
dge$NEURON_DIFFERENTIATION1 <- scale_zscore(dge$NEURON_DIFFERENTIATION1)
dge$NEUROPEPTIDE_RECEPTOR1 <- scale_zscore(dge$NEUROPEPTIDE_RECEPTOR1)
dge$EMT1 <- scale_zscore(dge$EMT1)

library(mgcv)
library(scales)
library(ggplot2)

### Define the GetGamTrends function for piecewise splines
GetGamTrends <- function(pdt_cells, genes, data, n_splines = 8) {
  win <- 30
  min_val <- 1e9
  gamx <- gamy <- gam_ci <- list()
  for (g in genes) {
    X <- as.matrix(pdt_cells)
    y <- as.vector(data[ ,g])
    # y <- data$ASCL1_Expression
    gam <- gam(y ~ s(X, bs = "cr", k = n_splines))
    XX <- seq(min(X) - 0.10, max(X) + 0.10, length.out = 500)
    YY <- predict(gam, newdata = data.frame(X = XX), type = "response")
    gamx[[g]] <- XX
    gamy[[g]] <- YY
    if (min(YY) < min_val) min_val <- min(YY)
    gam_ci[[g]] <- predict(gam, newdata = data.frame(X = XX), interval = "confidence", level = 0.95)
  }
  return(list(gamx = gamx, gamy = gamy, gam_ci = gam_ci))
}

### Filter the TP53/RB1 wt SCLC cells and SCLC-P cells
filter_genes <- grepl("^MT-|^MTMR|^MTND|NEAT1|TMSB4X|TMSB10|^RPS|^RPL|^RP11|^MRP|^FAU$|UBA52|MALAT|^IGH|^IGK|^IGL[CV]|^HBA|^HBB", rownames(dge))
filter_cells <- colnames(dge)[dge$SCLC_subtype_plus_TP53_RB1_wt %in% c('TP53/RB1-wt SCLC', 'SCLC-P')]

### Make the dataframe for plotting
Meta <- dge@meta.data
pval_df <- Meta[, c('pval_SCLC-A', 'pval_SCLC-N')]
pval_df <- pval_df / rowSums(pval_df)
pval_df <- pval_df[!(rownames(pval_df) %in% filter_cells),]

ind <- order(pval_df[, 'pval_SCLC-A'])

scaled01_imp_df <- Meta[,c("AXONOGENESIS1","NEURON_DIFFERENTIATION1",
                           "NEUROPEPTIDE_RECEPTOR1", "EMT1",
                           "ASCL1_Expression", "NEUROD_genes")]
scaled01_imp_df <- scaled01_imp_df[!(rownames(scaled01_imp_df) %in% filter_cells),]

pval_SCLCtype <- pval_df[ind, 'pval_SCLC-A']
scaled01_imp_df <- scaled01_imp_df[ind, ]

pal_GR <- colorRampPalette(c("red", "green"))
bin_SCLCtype <- cut(pval_SCLCtype, breaks = 500, labels = FALSE)
lut <- pal_GR(500)[bin_SCLCtype]
row_colors2 <- lut

### Plotting the first two
pdt <- pval_df[ind, 'pval_SCLC-A']
pdt2 <- as.integer(rank(pdt, ties.method = "min"))
gam_trends <- GetGamTrends(pdt2, c('ASCL1_Expression', 'NEUROD_genes'), scaled01_imp_df)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1))
for (i in 1:2) {
  g <- c('ASCL1_Expression', 'NEUROD_genes')[i]
  col <- c('red', 'green')[i]
  plot(pdt2[rev(seq_along(pdt2))], scaled01_imp_df[, g], pch = 20, col = 'lightgray', cex = 0.5,
       xlim = range(pdt2), ylim = c(-0.05, 1.05), xlab = 'Pvalue_A_vs_N', ylab = g,
       main = g)
  lines(gam_trends$gamx[[g]][rev(seq_along(gam_trends$gamx[[g]]))], gam_trends$gamy[[g]],
        col = col, lwd = 2)
  # legend("right", legend = "Binary Classification", title = "SCLC Type",
  #        fill = pal_GR(500), bty = "n", ncol = 2, cex = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


### Plotting the middle two
pdt <- pval_df[ind, 'pval_SCLC-A']
pdt2 <- as.integer(rank(pdt, ties.method = "min"))
gam_trends <- GetGamTrends(pdt2, c('AXONOGENESIS1', 'NEURON_DIFFERENTIATION1'), scaled01_imp_df)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1))
for (i in 1:2) {
  g <- c('AXONOGENESIS1', 'NEURON_DIFFERENTIATION1')[i]
  col <- c('red', 'green')[i]
  plot(pdt2[rev(seq_along(pdt2))], scaled01_imp_df[, g], pch = 20, col = 'lightgray', cex = 0.5,
       xlim = range(pdt2), ylim = c(-0.05, 1.05), xlab = 'Pvalue_A_vs_N', ylab = g,
       main = g)
  lines(gam_trends$gamx[[g]][rev(seq_along(gam_trends$gamx[[g]]))], gam_trends$gamy[[g]],
        col = col, lwd = 2)
  # legend("right", legend = "Binary Classification", title = "SCLC Type",
  #        fill = pal_GR(500), bty = "n", ncol = 2, cex = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

### Plotting the bottom two
pdt <- pval_df[ind, 'pval_SCLC-A']
pdt2 <- as.integer(rank(pdt, ties.method = "min"))
gam_trends <- GetGamTrends(pdt2, c('NEUROPEPTIDE_RECEPTOR1', 'EMT1'), scaled01_imp_df)
par(mfrow = c(1, 2), mar = c(5, 4, 2, 1))
for (i in 1:2) {
  g <- c('NEUROPEPTIDE_RECEPTOR1', 'EMT1')[i]
  col <- c('red', 'green')[i]
  plot(pdt2[rev(seq_along(pdt2))], scaled01_imp_df[, g], pch = 20, col = 'lightgray', cex = 0.5,
       xlim = range(pdt2), ylim = c(-0.05, 1.05), xlab = 'Pvalue_A_vs_N', ylab = g,
       main = g)
  lines(gam_trends$gamx[[g]][rev(seq_along(gam_trends$gamx[[g]]))], gam_trends$gamy[[g]],
        col = col, lwd = 2)
  # legend("right", legend = "Binary Classification", title = "SCLC Type",
  #        fill = pal_GR(500), bty = "n", ncol = 2, cex = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


