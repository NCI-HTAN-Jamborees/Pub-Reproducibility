
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
library(jjb)



####Make Epithelial subset
options(future.globals.maxSize = 80000 * 1024^2)
setwd("/Volumes/Backup_Plus/HTAN/Epithelial/")
dge<-readRDS("../Overall_HTAN.rds")
dge<-subset(dge,cells=colnames(dge)[dge$cell_type_coarse=="Epithelial"])
library(stringr)
library("ggplot2")
library("dplyr")
library('reshape2')
#library(entropy)
library(tidyr)
library(plyr)
library(entropy)
library(ggrepel)
library(RColorBrewer)
library(pals)

Train_data<-read.table("./Epithelial_metadata.csv",sep=",",header = T,row.names = 1)
dge@meta.data<-cbind(dge@meta.data,Train_data[,43:46])


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

####Figure S2A
###Subtype Diversity and Frequency
subtype_freq = metadata_df %>% dplyr::group_by(donor_id, SCLC_subtype) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq = n/sum(n)) %>% 
  dplyr::select(-n) %>% pivot_wider(names_from = SCLC_subtype, values_from = freq) 
subtype_freq = as.data.frame(subtype_freq)
subtype_freq[is.na(subtype_freq)] <- 0
rownames(subtype_freq) = subtype_freq$donor_id
donor_id = subtype_freq$donor_id
subtype_freq = subtype_freq %>% dplyr::select(-donor_id)

SCLCmajor = colnames(subtype_freq)[apply(subtype_freq,1,which.max)]
names(SCLCmajor) = rownames(subtype_freq)

subtype_diversity = sort(apply(as.matrix(subtype_freq), 1, entropy))

plot_df0 = metadata_df %>% dplyr::select(c('donor_id','subtype_uncertainty'))
colnames(plot_df0) = c('Sample','Uncertainty')
plot_df0$SCLCtype_major = SCLCmajor[plot_df0$Sample]
plot_df0$Sample = factor(plot_df0$Sample, levels = names(subtype_diversity))

p0 = plot_df0%>% ggplot(aes(Sample,Uncertainty,fill=SCLCtype_major)) + geom_boxplot(outlier.size=-1) + xlab('SCLC subtype') + ylab('Subtype uncertainty') + theme_bw() + 
  theme(text = element_text(size=15), strip.text = element_text(size=25),axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill=guide_legend(title='Major\nSubclonal\nSubtype'))

tally_df = metadata_df %>% dplyr::select(donor_id,SCLC_subtype)
tally_df$donor_id = as.factor(tally_df$donor_id)
tally_df$SCLC_subtype = as.factor(tally_df$SCLC_subtype)
counts = tally_df %>% dplyr::group_by(donor_id, SCLC_subtype, .drop=F) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(freq = n/sum(n)) %>% 
  dplyr::select(-n) %>% pivot_wider(names_from = SCLC_subtype, values_from = freq) %>% dplyr::rename(Sample = donor_id)

plot_df1 = counts %>% melt(id = 'Sample') %>% mutate(Population_Frequency = 'Frequency')
plot_df1$SCLCmajor = SCLCmajor[as.character(plot_df1$Sample)]
plot_df1$Sample = factor(plot_df1$Sample, levels = names(subtype_diversity))
plot_df1$SCLCmajor = factor(plot_df1$SCLCmajor, levels = sort(unique(as.character(plot_df1$SCLCmajor))) )
p1 = plot_df1 %>% ggplot(aes(Sample,value,fill=variable)) + geom_bar(position='fill', stat='identity') + 
  xlab('Sample') + ylab('Population Fraction') + theme_bw() + theme(text = element_text(size=15), axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(fill = 'Subtype')

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(p0), ggplotGrob(p1), size = "last"))

####Done
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


n_pc=100

dge <- RunUMAP(dge, dims = 1:100)


DimPlot(dge, reduction = "umap",group.by="cell_type_fine",label=T,raster = T)
ggsave(file="Umap_fineID.pdf",width = 20,height = 15)

DimPlot(dge, reduction = "umap",group.by="donor_id",label=T,raster = T)
ggsave(file="Umap_donor.pdf",width = 20,height = 15)

####Figure 1H
p1<-FeaturePlot(dge,features=c("ASCL1"),cols = c("purple","red","yellow"))
p2<-FeaturePlot(dge,features=c("NEUROD1"),cols = c("purple","red","yellow"))
p3<-FeaturePlot(dge,features=c("POU2F3"),cols = c("purple","red","yellow"))
p4<-FeaturePlot(dge,features=c("YAP1"),cols = c("purple","red","yellow"))
grid.arrange(p1, p2, p3, p4, ncol = 2)

###Make the tumor object
dge<-subset(Epithelial,cells=colnames(Epithelial)[Epithelial$cell_type_fine %in% c("SCLC-A","SCLC-N","SCLC-P")])

####Figure 2A
dge<-SetIdent(dge,value = as.vector(dge$cell_type_fine))
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


