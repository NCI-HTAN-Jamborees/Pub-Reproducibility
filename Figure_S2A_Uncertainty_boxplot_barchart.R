
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

####Figure S2A
###Subtype Diversity and Frequency
dge<-readRDS("./SCLC_Epithelial.rds")

metadata_df = dge@meta.data
metadata_df$batch = as.character(metadata_df$batch)
metadata_df$batch[grepl('RU1108',metadata_df$batch)] = 'RU1108a'
#### Notice: might need to change "." to "-"
metadata_df$subtype_uncertainty = apply(metadata_df[,c("pval_SCLC.A","pval_SCLC.N","pval_SCLC.P")], 1, function(x) entropy(x))
any(is.nan(metadata_df$subtype_uncertainty))
                                        
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

