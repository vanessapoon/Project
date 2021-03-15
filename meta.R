setwd("~/OneDrive - University of Birmingham/Year 3/Research Project/data") 

library(tidyverse)
# install.packages("cowplot")
# install.packages("hrbrthemes")
library(cowplot)
library(hrbrthemes)

results_all <- readRDS("results_all.RDS")

names(results_all) <- plyr::mapvalues(x = names(results_all), from = "GSE146028_IFNg", to = "GSE146028_IFNg+LPS")

datasets <- read.delim("Datasets.txt")

sub <- datasets$ID [c(1,3,5,8)]

lfc <- list()
for (i in sub){
  lfc[[i]] <- results_all[[i]]$log2FoldChange
}

padj <- list()
for (i in sub){
  padj[[i]] <- results_all[[i]]$padj
}

genes <- rownames(results_all[[1]])

lfc_df <- lfc %>% as.data.frame()
rownames(lfc_df) <- genes

padj_df <- padj %>% as.data.frame()
rownames(padj_df) <- genes

library(org.Hs.eg.db)
Symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = substr(genes, 1, 15), keytype = "ENSEMBL", column = "SYMBOL")

targets <- c("CHRD", "CHRDL1", "DAND5", "FST", "FSTL1", "FSTL3", "GREM1", "GREM2", "NBL1", "NOG", "SOST", "SOSTDC1", "TWSG1")
# targets <- c("BMP2", "BMP4", "BMP7", "BMPR1A", "BMPR1B", "BMPR2", "BMPER", "ENG", "MIF", "CD44", "CD74", "MIF", "KDR")

# heatmap for BMP antagonists 
lfc_select <- lfc_df [Symbol %in% targets, ]
rownames(lfc_select) <- Symbol [Symbol %in% targets]
lfc_long <- lfc_select %>% mutate(gene_id = factor(rownames(lfc_select))) %>% pivot_longer(cols = 1:3)
padj_select <- padj_df [Symbol %in% targets, ]
rownames(padj_select) <- Symbol [Symbol %in% targets]
padj_long <- padj_select %>% mutate(gene_id = rownames(padj_select)) %>% pivot_longer(cols = 1:3)
padj_long$value [is.na(padj_long$value)] <- 1

lfc_long$padj <- padj_long$value

lfc_long$signif <- case_when(
  lfc_long$padj < 0.001 ~ "***",
  lfc_long$padj < 0.01 ~ "**",
  lfc_long$padj < 0.05 ~ "*") 

max <- max(abs(lfc_long$value), na.rm = T)
p1 <- lfc_long %>% 
  ggplot(aes(x = name, y = gene_id, fill = value))+
  geom_tile()+
  geom_text(aes(label = signif))+
  labs(title = "Differential gene expression of BMP antagonists in LPS stimulated macrophages", x ="", y = "", 
       fill = expression(paste("log" [2], "-fc")))+
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-max, max)) +
  scale_y_discrete(limits = rev(levels(lfc_long$gene_id)))+
  theme_minimal()+
  theme(line = element_line(colour = "black", size = 0.5),
        text = element_text(family = "serif"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.background = element_blank(), panel.grid = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size=10, face = "bold", angle = 90, hjust = 1),
        axis.text.y = element_text(size=10, face = "bold"))
p1

# pdf(file = "DEG analysis for BMP antagonists in HSC.pdf", width=5, height=5)
# p1
# dev.off()

# check for correlation of datasets ####
lfc_mat <- as.matrix(lfc_select)
corrmat <- Hmisc::rcorr(lfc_mat)

install.packages("corrplot"")
corrplot::corrplot(corrmat$r)

pheatmap::pheatmap(corrmat$r)

