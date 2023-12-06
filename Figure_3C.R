###To use the published cellphoneDB result
###https://github.com/dpeerlab/SCLC_atlas-HTAN/raw/main/out.SCLC.060121/cellphonedb.DE.homotypic_tumor.SCLC-A_vs_SCLC-N.txt


###Otherwise, start with the SCLC object and make the input for cellphoneDB
###Load the SCLC Epithelial cell object
###This dataset can be found here: https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec

####Figure 2C CellphoneDB result
library(Seurat)
dge <- readRDS("./Download_SCLC_Epithelial.rds")

####Change geneID to gene symbols
Features <- dge@assays$RNA@meta.features
RNA <- dge@assays$RNA
RNA@counts@Dimnames[[1]] <- as.vector(Features$feature_name)
RNA@data@Dimnames[[1]] <- as.vector(Features$feature_name)
dge@assays$RNA <- RNA

dge <- SetIdent(dge, value = as.vector(dge$SCLC_subtype))

####For computational efficiency
dge <- subset(dge, idents = c("SCLC-N","SCLC-A"), downsample = 500)

####For CellphoneDB, needs to re-normalize using RC instead of log normalization. 
dge <- NormalizeData(object = dge, normalization.method = "RC", scale.factor = 1e6)
dge_mtx <- dge@assays$RNA@data
write.table(dge_mtx, "./SCLC_A_N_count.txt", sep = "\t", row.names = T, col.names = T, quote = F)
DF <- NULL
DF$Cell <- colnames(dge)
DF$cell_type <- as.vector(dge@active.ident)
write.table(DF, "./SCLC_A_N_meta.txt", sep = "\t", row.names = F, col.names = T, quote = F)

####Next, add "Gene" in the count file. 
####Create a folder called SCLC_A_N for cellphonedb output
####Run this command for cellphoneDB analysis:
####cellphonedb method statistical_analysis ./SCLC_A_N_meta.txt ./SCLC_A_N_count.txt --output-path=./SCLC_A_N --counts-data=hgnc_symbol

####Finish the CPDB analysis, Now processing the result for visualization in dotplot
DF <- read.table("./SCLC_A_N/pvalues.txt", sep = "\t", header = T)
####Since only self-interaction is wanted, only subset those pairs
Sig_columnpair <- c("SCLC.A.SCLC.A", "SCLC.N.SCLC.N")
write.table(Sig_columnpair, "./SCLC_A_N/columns_direction.txt", sep = "\n",col.names = F,row.names = F,quote = F)
rownames(DF) <- make.unique(as.vector(DF$interacting_pair))
DF_orig <- DF
DF <- DF[,colnames(DF) %in% Sig_columnpair]
Sig_rows <- NULL
for (i in 1:length(Sig_columnpair)) {
  Sig_rows_temp <- rownames(DF)[DF[,i]<0.05]
  if (length(Sig_rows_temp) >= 1) {
    Sig_rows <- union(Sig_rows,Sig_rows_temp)}
}
write.table(Sig_rows, "./SCLC_A_N/rows_direction.txt", sep = "\n",col.names = F,row.names = F,quote = F)

####Notice: Go in the columns_direction.txt file and change the . to | which is required by CPDB
# cellphonedb plot dot_plot --output-path=./SCLC_A_N --pvalues-path=./SCLC_A_N/pvalues.txt --means-path=./SCLC_A_N/means.txt --rows ./SCLC_A_N/rows_direction.txt --columns ./SCLC_A_N/columns_direction.txt
####Subset the row direction file with the interactions in the paper

### An attempt to make the volcano plot in the manuscript but failed 
# ####Visualize in a volcano-plot type of figure
# pvalue <- read.table("./SCLC_A_N/pvalues.txt", sep = "\t", header = T)
# means <- read.table("./SCLC_A_N/means.txt", sep = "\t", header = T)
# DF_Volcano <- cbind(pvalue$interacting_pair, pvalue$SCLC.A.SCLC.A, pvalue$SCLC.N.SCLC.N, 
#                     means$SCLC.A.SCLC.A, means$SCLC.N.SCLC.N)
# colnames(DF_Volcano) <- c("Pairs", "Pvalue_SCLC_A", "Pvalue_SCLC_N", "Expression_SCLC_A", "Expression_SCLC_N")
# 
# DF_Volcano$Pvalue_SCLC_A <- as.numeric(DF_Volcano$Pvalue_SCLC_A)
# DF_Volcano$Pvalue_SCLC_N <- as.numeric(DF_Volcano$Pvalue_SCLC_N)
# DF_Volcano$Expression_SCLC_A <- as.numeric(DF_Volcano$Expression_SCLC_A)
# DF_Volcano$Expression_SCLC_N <- as.numeric(DF_Volcano$Expression_SCLC_N)
# DF_Volcano<-as.data.frame(DF_Volcano)
# 
# top_pairs <- head(DF_Volcano[order(DF_Volcano$Pvalue_SCLC_N), "Pairs"], 20)
# 
# # Create a volcano plot
# volcano_plot <- ggplot(DF_Volcano, aes(x = Expression_SCLC_N - Expression_SCLC_A, y = -log10(Pvalue_SCLC_N))) +
#   geom_point(aes(color = factor(sign(Expression_SCLC_N - Expression_SCLC_A))),
#              size = 2, alpha = 0.8) +
#   #scale_color_manual(values = c("red", "blue")) +  # Specify colors for upregulated and downregulated
#   labs(x = "Expression Difference (Right: SCLC_N, Left: SCLC_A)", y = "-log10(P-value)") +
#   theme_minimal() +
#   geom_text_repel(data = subset(DF_Volcano, Pairs %in% top_pairs),
#                   aes(label = Pairs), box.padding = 0.5, size = 3)+NoLegend()  # Add labels for the top 20 pairs with repulsion
# 
# # Show the plot
# print(volcano_plot)


