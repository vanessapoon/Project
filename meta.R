setwd("~/OneDrive - University of Birmingham/Year 3/Research Project/data") 

library(tidyverse)
library(cowplot)
library(hrbrthemes)

results_all <- readRDS("results_all.RDS")

names(results_all) <- plyr::mapvalues(x = names(results_all), from = "GSE146028_IFNg", to = "GSE146028_IFNg+LPS")

datasets <- read.delim("Datasets.txt")

sub <- datasets$ID

#creates new lists of just log2FC and padj from results_all
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

# selected gene targets
# classic m1/ m2 markers
targets  <- c("TNF", "IL6", "NOS2", "IL12B", "IL23A", "CXCL9", "CXCL10", "CXCL11", "IL10", "ARG1", "CCL17", "CCL22", "TGFB1", "CD14", "CD16", "NOS2", "CD163", "CD80", "MRC1", "ITGAM", "IRF5", "IFNG", "STAT1", "IRF4", "STAT3") 
# bmp antagonists 
# targets <- c("CER1","CHRD", "CHRDL1", "DAND5", "FST", "FSTL1", "FSTL3", "GREM1", "GREM2", "NBL1", "NOG", "SOST", "SOSTDC1", "TWSG1")
# bmp pathway
# targets <- c("BAMBI", "BMP2", "BMP4", "BMP6", "BMP7", "BMPR1A", "BMPR1B", "BMPR2", "MAP3K","MAPK1", "SMAD1", "SMAD4","SMAD5", "SMAD6", "SMAD7" , "SMAD9")
# TGFB pathway
# targets <- c("APC", "CDH1", "CREBBP", "EP300", "MAP2K1", "MAP3K7", "MAPK3", "SKIL", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "TAB1", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "ZFYVE9")
# activin receptor pathway
# targets  <- c("ACVR1", "ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "ACVRL1", "BMP10", "BMPR2", "CER1", "CFC1", "CFC1B", "SMAD2", "SMAD3", "SMAD7", "DAND5", "FGF10", "GDF2", "GDF6", "GDF7", "NODAL")


#heatmap for selected genes
#column number needs to match number of datsets for pivot_longer
lfc_select <- lfc_df [Symbol %in% targets, ]
rownames(lfc_select) <- Symbol [Symbol %in% targets]
lfc_long <- lfc_select %>% mutate(gene_id = factor(rownames(lfc_select))) %>% pivot_longer(cols = 1:8)
padj_select <- padj_df [Symbol %in% targets, ]
rownames(padj_select) <- Symbol [Symbol %in% targets]
padj_long <- padj_select %>% mutate(gene_id = rownames(padj_select)) %>% pivot_longer(cols = 1:8)
padj_long$value [is.na(padj_long$value)] <- 1

lfc_long$padj <- padj_long$value

lfc_long$signif <- case_when(
  lfc_long$padj < 0.001 ~ "***",
  lfc_long$padj < 0.01 ~ "**",
  lfc_long$padj < 0.05 ~ "*") 

max <- max(abs(lfc_long$value), na.rm = T)
p1 <- lfc_long %>% 
  ggplot(aes(x = factor(name, level = c("GSE146028_IFNg.LPS", "GSE98368_IFNg", "GSE100382_LPS", "GSE163165_LPS", "GSE135753_LPS", "GSE146028_IL4", "GSE132732_IL4", "GSE123603_IL4")), y = gene_id, fill = value))+
  geom_tile()+
  geom_text(aes(label = signif))+
  labs(title = "", x ="", y = "", 
       fill = expression(paste("log" [2], "-fc")))+
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-max, max)) +
  scale_y_discrete(limits = rev(levels(lfc_long$gene_id)))+
  theme_minimal()+
  theme(line = element_line(colour = "black", size = 0.5),
        #text = element_text(family = "serif"),
        # plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.background = element_blank(), panel.grid = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.title = element_text(size = 10, hjust = 0.5),
        axis.text.x = element_text(size=10, angle = 90, hjust = 1),
        axis.text.y = element_text(size=10))
p1

#check for correlation of datasets
lfc_mat <- as.matrix(lfc_select)
corrmat <- Hmisc::rcorr(lfc_mat)

corrplot::corrplot(corrmat$r)

pheatmap::pheatmap(corrmat$r)


# META FOR LPS
library(metaRNASeq)
subLPS <- datasets$ID[datasets$Treatment == "LPS"]

results_meta <- list()

for(i in subLPS){  
  results_meta[[i]] <- results_all[[i]]
}

rawpval <- list()
FC <- list()
adjpval <- list()

for(i in names(results_meta)){
  rawpval[[i]] <- results_meta[[i]]$pvalue
  FC[[i]] <- results_meta[[i]]$log2FoldChange
  adjpval[[i]] <- results_meta[[i]]$padj
}

DE_LPS <- mapply(adjpval, FUN = function(x) ifelse (x <= 0.05, 1, 0))

filtered <- lapply(adjpval, FUN=function(pval) which(is.na(pval)))

for(i in names(rawpval)){
  rawpval[[i]][filtered[[i]]] <- NA
}
#correlation of datasets
FC_mat <- data.frame(FC) %>% as.matrix()

corrmat <- Hmisc::rcorr(FC_mat)
pheatmap::pheatmap(corrmat$r)

fishcomb_LPS <- fishercomb(rawpval, BHth = 0.05)
sum(fishcomb_LPS$adjpval < 0.01, na.rm = T)

#need to specify number of replicates for each study
#LPS = 2,3,3 IL4 = 2,5,3
invnormcomb_LPS <- invnorm(rawpval,nrep=c(2,3,3), BHth = 0.05)  #save separately
sum(invnormcomb_LPS$adjpval < 0.01, na.rm = T) #save separately

DEresults_LPS <- data.frame(DE_LPS,
                            "DE.fishercomb"=ifelse(fishcomb_LPS$adjpval<=0.01,1,0),
                            "DE.invnorm"=ifelse(invnormcomb_LPS$adjpval<=0.01,1,0)) #rename DEresults object
head(DEresults_LPS)

#remove contradictory differential expression
signsFC <- mapply(FC, FUN=function(x) sign(x))
sumsigns <- apply(signsFC, 1, sum)
commonsgnFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns), 0)

signsFC <- mapply(FC, FUN=function(x) sign(x))
signsFC[DE_LPS == 0] <- 0
sumsigns <- apply(signsFC, 1, sum, na.rm = T)
commonsgnFC <- ifelse(abs(sumsigns)==rowSums(DE_LPS, na.rm = T), sign(sumsigns), 0)

sum(commonsgnFC !=0, na.rm =T)

DEresults_LPS$commonsgn <- commonsgnFC

# extract ENSEMBL IDs for DEGs, either all or only up/down-regulated genes
DE_genes_all_LPS <- names(Symbol) [DEresults_LPS$DE.fishercomb == 1 & DEresults_LPS$DE.invnorm == 1 & DEresults_LPS$commonsgn!= 0 & !is.na(DEresults_LPS$commonsgn)]
DE_genes_up_LPS <- names(Symbol) [DEresults_LPS$DE.fishercomb == 1 & DEresults_LPS$DE.invnorm == 1 & DEresults_LPS$commonsgn == 1 & !is.na(DEresults_LPS$commonsgn)]
DE_genes_down_LPS <- names(Symbol) [DEresults_LPS$DE.fishercomb == 1 & DEresults_LPS$DE.invnorm == 1 & DEresults_LPS$commonsgn == -1 & !is.na(DEresults_LPS$commonsgn)]

# META FOR IL4
library(metaRNASeq)
subIL4 <- datasets$ID[datasets$Treatment == "IL4"]
results_meta <- list()

for(i in subIL4){  
  results_meta[[i]] <- results_all[[i]]
}

rawpval <- list()
FC <- list()
adjpval <- list()

for(i in names(results_meta)){
  rawpval[[i]] <- results_meta[[i]]$pvalue
  FC[[i]] <- results_meta[[i]]$log2FoldChange
  adjpval[[i]] <- results_meta[[i]]$padj
}

DE_IL4 <- mapply(adjpval, FUN = function(x) ifelse (x <= 0.05, 1, 0))

filtered <- lapply(adjpval, FUN=function(pval) which(is.na(pval)))

for(i in names(rawpval)){
  rawpval[[i]][filtered[[i]]] <- NA
}
#correlation of datasets
FC_mat <- data.frame(FC) %>% as.matrix()

corrmat <- Hmisc::rcorr(FC_mat)
pheatmap::pheatmap(corrmat$r)

fishcomb_IL4 <- fishercomb(rawpval, BHth = 0.05)
sum(fishcomb_IL4$adjpval < 0.01, na.rm = T)

#need to specify number of replicates for each study
#LPS = 2,3,3 IL4 = 2,5,3
invnormcomb_IL4 <- invnorm(rawpval,nrep=c(2,5,3), BHth = 0.05)  #save separately
sum(invnormcomb_IL4$adjpval < 0.01, na.rm = T) #save separately

DEresults_IL4 <- data.frame(DE_IL4,
                            "DE.fishercomb"=ifelse(fishcomb_IL4$adjpval<=0.01,1,0),
                            "DE.invnorm"=ifelse(invnormcomb_IL4$adjpval<=0.01,1,0)) #rename DEresults object
head(DEresults_IL4)

#remove contradictory differential expression
signsFC <- mapply(FC, FUN=function(x) sign(x))
sumsigns <- apply(signsFC, 1, sum)
commonsgnFC <- ifelse(abs(sumsigns)==dim(signsFC)[2], sign(sumsigns), 0)

signsFC <- mapply(FC, FUN=function(x) sign(x))
signsFC[DE_IL4 == 0] <- 0
sumsigns <- apply(signsFC, 1, sum, na.rm = T)
commonsgnFC <- ifelse(abs(sumsigns)==rowSums(DE_IL4, na.rm = T), sign(sumsigns), 0)

sum(commonsgnFC !=0, na.rm =T)

DEresults_IL4$commonsgn <- commonsgnFC

# extract ENSEMBL IDs for DEGs, either all or only up/down-regulated genes
DE_genes_all_IL4 <- names(Symbol) [DEresults_IL4$DE.fishercomb == 1 & DEresults_IL4$DE.invnorm == 1 & DEresults_IL4$commonsgn!= 0 & !is.na(DEresults_IL4$commonsgn)]
DE_genes_up_IL4 <- names(Symbol) [DEresults_IL4$DE.fishercomb == 1 & DEresults_IL4$DE.invnorm == 1 & DEresults_IL4$commonsgn == 1 & !is.na(DEresults_IL4$commonsgn)]
DE_genes_down_IL4 <- names(Symbol) [DEresults_IL4$DE.fishercomb == 1 & DEresults_IL4$DE.invnorm == 1 & DEresults_IL4$commonsgn == -1 & !is.na(DEresults_IL4$commonsgn)]

# meta and other 2 data sets altogether

IFNg <- results_all$GSE98368_IFNg
IFNgandLPS <- results_all$`GSE146028_IFNg+LPS`

IFNg <-  IFNg %>% dplyr::select(log2FoldChange, padj) 

#subsetting up and downregulated genes
IFNg_up<- dplyr::filter(IFNg, log2FoldChange > 0, !is.na(padj))
IFNg_down <-  dplyr::filter(IFNg, log2FoldChange < 0, !is.na(padj))

IFNgandLPS <- IFNgandLPS %>% dplyr::select(log2FoldChange, padj)

IFNgandLPS_up <-  dplyr::filter(IFNgandLPS, log2FoldChange > 0, !is.na(padj))
IFNgandLPS_down <-  dplyr::filter(IFNgandLPS, log2FoldChange < 0, !is.na(padj))

library(AnnotationDbi)
library(org.Hs.eg.db)
library(GenomicFeatures)

# gtf <- "gencode.gtf"
# txdb <- makeTxDbFromGFF(gtf)
# saveDb(txdb,"txdb.filename")

# txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
# tab <- table(txdf$GENEID)
# txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

#adding gene symbols to gene lists
IFNg_up <-as.data.frame(IFNg_up) %>%
  mutate(genes = rownames(IFNg_up)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

IFNg_down <- as.data.frame(IFNg_down) %>%
  mutate(genes = rownames(IFNg_down)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

IFNgandLPS_up <- as.data.frame(IFNgandLPS_up) %>%
  mutate(genes = rownames(IFNgandLPS_up)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

IFNgandLPS_down <- as.data.frame(IFNgandLPS_down) %>%
  mutate(genes = rownames(IFNgandLPS_down)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

IFNg_up_gene_list <- IFNg_up$genes [IFNg_up$padj < 0.01] 
IFNg_up_gene_list <- substr(IFNg_up_gene_list, 1, 15) 
IFNg_up_gene_list

IFNg_down_gene_list <- IFNg_down$genes [IFNg_down$padj < 0.01] 
IFNg_down_gene_list <- substr(IFNg_down_gene_list, 1, 15) 
IFNg_down_gene_list

IFNgandLPS_up_gene_list <- IFNgandLPS_up$genes [IFNgandLPS_up$padj < 0.01] 
IFNgandLPS_up_gene_list <- substr(IFNgandLPS_up_gene_list, 1, 15) 
IFNgandLPS_up_gene_list

IFNgandLPS_down_gene_list <- IFNgandLPS_down$genes [IFNgandLPS_down$padj < 0.01] 
IFNgandLPS_down_gene_list <- substr(IFNgandLPS_down_gene_list, 1, 15) 
IFNgandLPS_down_gene_list

library(clusterProfiler)
enr_GO <- enrichGO(DE_genes_up_IL4, org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
enr_GO <- setReadable(enr_GO, org.Hs.eg.db, keyType = "ENSEMBL") #converting to gene symbols
head(enr_GO)

#barplot to show enriched pathways
library(enrichplot)
barplot(enr_GO, showCategory=10, label_format=21) +
  labs(title = "IL4", x = "Counts") + 
  expand_limits(x = c(0, 150)) +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.text.y = element_text(size=10)) 



#adding gene symbols 
LPS_invnorm <- as.data.frame(invnormcomb_LPS$adjpval) %>%
  mutate(genes = rownames(DE_genes_all_LPS)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

LPS_fish <- as.data.frame(fishcomb_LPS$adjpval) %>%
  mutate(genes = rownames(DE_genes_all_LPS)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

IL4_invnorm <- as.data.frame(invnormcomb_IL4$adjpval) %>%
  mutate(genes = rownames(DE_genes_all_IL4)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

IL4_fish <- as.data.frame(fishcomb_IL4$adjpval) %>%
  mutate(genes = rownames(DE_genes_all_IL4)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(genes, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"))

# selecting minimum p values for all genes to show on heatmap 

# DGE <- list()
# saveRDS(DGE, file = "DGE.RDS")

DGE <- readRDS("DGE.RDS")
DGE[["IL4_fishcomb_padj"]] <-  fishcomb_IL4$adjpval
saveRDS(DGE, file = "DGE.RDS")

DGE <- readRDS("DGE.RDS")
DGE[["IL4_invnorm_padj"]] <-  invnormcomb_IL4$adjpval
saveRDS(DGE, file = "DGE.RDS")

DGE <- readRDS("DGE.RDS")
DGE[["LPS_fishcomb_padj"]] <-  fishcomb_LPS$adjpval
saveRDS(DGE, file = "DGE.RDS")

DGE <- readRDS("DGE.RDS")
DGE[["LPS_invnorm_padj"]] <-  invnormcomb_LPS$adjpval
saveRDS(DGE, file = "DGE.RDS")

DGE <- readRDS("DGE.RDS")
DGE[["IFNg"]] <-  IFNg$padj
saveRDS(DGE, file = "DGE.RDS")

DGE <- readRDS("DGE.RDS")
DGE[["IFNgandLPS"]] <-  IFNgandLPS$padj
saveRDS(DGE, file = "DGE.RDS")

genes <- rownames(results_all[[1]])   #somehow getting the gene ids (load file first)

allpadj_df <- DGE %>% as.data.frame()
rownames(allpadj_df) <- genes
library(tidyverse)
library(org.Hs.eg.db)
Symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = substr(genes, 1, 15), keytype = "ENSEMBL", column = "SYMBOL")

minpadj_df <- as.data.frame(apply(allpadj_df, 1, min, na.rm=TRUE))

names(minpadj_df)[names(minpadj_df) == "apply(allpadj_df, 1, min, na.rm = TRUE)"] <- "minimum_padj"

minpadj_df <- minpadj_df %>% filter(minimum_padj > 0) %>% arrange(minimum_padj) %>% top_n(-2000) 

#creating matrix of log2foldchagnes of top 1000 regulated genes for creating overview heatmap
matrix_heatmap_overview <- lfc_df [rownames(lfc_df) %in% rownames(minpadj_df), ] %>% as.matrix()
library(pheatmap)
library(ComplexHeatmap)
library(dendsort)
library(dendextend)
library(circlize)
library(RColorBrewer)
library(viridis)
library(gridtext)

mat_cluster_rows <- hclust(dist(matrix_heatmap_overview), method = "ward.D")
# mat_cluster_rows <- color_branches(mat_cluster_rows, k = 4)
plot(mat_cluster_rows, main = "Unsorted Dendrogram", xlab = "", sub = "")

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...), isReverse=TRUE, type="min"))

plot(mat_cluster_rows, main = "Sorted Dendrogram", xlab = "", sub = "")
mat_cluster_cols <- sort_hclust(hclust(dist(t(matrix_heatmap_overview))), method = "ward.D")

row_cluster <- cutree(mat_cluster_rows, k = 4)
row_cluster <- plyr::mapvalues(row_cluster, from = c(1,2,3,4), to = c(2,4,1,3))
row_cluster <- factor(row_cluster)

# Data frame with row cluster annotations.
mat_row <- data.frame(cluster = as.factor(row_cluster))
rownames(mat_row) <- names(row_cluster)

# List with colors for each annotation.
mat_colors <- list(cluster = brewer.pal(length(unique(row_cluster)), "YlGn"))
names(mat_colors$cluster) <- unique(row_cluster)

# pheatmap::pheatmap(matrix_heatmap_overview,
# show_rownames = FALSE,
# cluster_cols = mat_cluster_cols,
# cluster_rows = mat_cluster_rows,
#  drop_levels = TRUE,
# annotation_row = mat_row,
#  annotation_color = mat_colors,
#   scale = "row")

enr_GO_clust <- list()
enr_GO_clust_df <- list()

for (i in 1:length(unique(row_cluster))){
  geneList <- substr(names(row_cluster) [row_cluster == i], 1, 15)
  enr_GO_clust[[i]] <- enrichGO(geneList, ont = "BP", OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                                pvalueCutoff = 0.05)
  enr_GO_clust_df[[i]] <- head(enr_GO_clust[[i]], n = "Inf")
}

enr_KEGG_clust <- list()
enr_KEGG_clust_df <- list()

for (i in 1:length(unique(row_cluster))){
  geneList <- substr(names(row_cluster) [row_cluster == i], 1, 15)
  geneList <- mapIds(org.Hs.eg.db, geneList, keytype = "ENSEMBL", column = "ENTREZID",
                     pvalueCutoff = 0.05)
  enr_KEGG_clust[[i]] <- enrichKEGG(geneList, organism = "hsa")
  enr_KEGG_clust_df[[i]] <- head(enr_KEGG_clust[[i]], n = "Inf")
}

text_list = list(
  text1 = "complement activation \nresponse to IFNg & LPS \nmacrophage activation",
  text2 = "vesicle organisation \nRas and Rab signal transduction \nregulation of cell size",
  text3 = "autophagy",
  text4 = "metabolic regulation \nleukocyte differentiation \nECM organisation \nERK1 and ERK2 signaling"
)

ha = HeatmapAnnotation(bar = anno_empty(row_cluster, border = FALSE,
                                        width = max_text_width(unlist(text_list)) + unit(4, "mm")),
                       which = "row")

Heatmap(matrix_heatmap_overview,
        cluster_rows = mat_cluster_rows,
        cluster_columns = mat_cluster_cols,
        heatmap_legend_param = list(title = expression(paste("log" [2], "-fc"))),
        right_annotation = ha,
        row_split = 4,
        show_row_names = FALSE)
for(i in 1:4) {
  decorate_annotation("bar", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
  })
}

#heatmap for selected genes with ComplexHeatmap package####


# selected gene targets
# classic m1/ m2 markers
targets  <- c("TNF", "IL6", "NOS2", "IL12B", "IL23A", "CXCL9", "CXCL10", "CXCL11", "IL10", "ARG1", "CCL17", "CCL22", "TGFB1", "CD14", "CD16", "NOS2", "CD163", "CD80", "MRC1", "ITGAM", "IRF5", "IFNG", "STAT1", "IRF4", "STAT3") 
# bmp antagonists 
# targets <- c("CER1","CHRD", "CHRDL1", "DAND5", "FST", "FSTL1", "FSTL3", "GREM1", "GREM2", "NBL1", "NOG", "SOST", "SOSTDC1", "TWSG1")
# bmp pathway
# targets <- c("BAMBI", "BMP2", "BMP4", "BMP6", "BMP7", "BMPR1A", "BMPR1B", "BMPR2", "MAP3K","MAPK1", "SMAD1", "SMAD4","SMAD5", "SMAD6", "SMAD7" , "SMAD9")
# TGFB pathway
# targets <- c("APC", "CDH1", "CREBBP", "EP300", "MAP2K1", "MAP3K7", "MAPK3", "SKIL", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "TAB1", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "ZFYVE9")
# activin receptor pathway
# targets  <- c("ACVR1", "ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "ACVRL1", "BMP10", "BMPR2", "CER1", "CFC1", "CFC1B", "SMAD2", "SMAD3", "SMAD7", "DAND5", "FGF10", "GDF2", "GDF6", "GDF7", "NODAL")
# targets <- c("GREM1")

lfc_select <- lfc_df [Symbol %in% targets, ] %>% as.matrix()
rownames(lfc_select) <- Symbol [Symbol %in% targets]
padj_select <- padj_df [Symbol %in% targets, ] %>% as.matrix()
rownames(padj_select) <- Symbol [Symbol %in% targets]


padj_asterisk <- matrix(case_when(
  padj_select < 0.001 ~ "***",
  padj_select < 0.01 ~ "**",
  padj_select < 0.05 ~ "*"), ncol = ncol(padj_select))

meta <- data.frame(LPS = case_when(apply(DEresults_LPS[,4:6], 1, sum) == 3 ~ apply(cbind(fishcomb_LPS$adjpval, invnormcomb_LPS$adjpval), 1, max),
                                   apply(DEresults_LPS[,4:6], 1, sum) != 3 ~ 1),
                   IL4 = case_when(apply(DEresults_IL4[,4:6], 1, sum) == 3 ~ apply(cbind(fishcomb_IL4$adjpval, invnormcomb_IL4$adjpval), 1, max),
                                   apply(DEresults_IL4[,4:6], 1, sum) != 3 ~ 1)
)
meta <- meta [Symbol %in% targets, ]
rownames(meta) <- Symbol [Symbol %in% targets]

ha = rowAnnotation("LPS" = anno_text(case_when(
  meta$LPS < 0.001 ~ "***",
  meta$LPS < 0.01 ~ "**",
  meta$LPS < 0.05 ~ "*",
  meta$LPS >= 0.05 ~ "",
  is.na(meta$LPS) ~ "" ),
  width = unit(6.5, "mm")),
  "IL4" = anno_text(case_when(
    meta$IL4 < 0.001 ~ "***",
    meta$IL4 < 0.01 ~ "**",
    meta$IL4 < 0.05 ~ "*",
    meta$IL4 >=0.05 ~ "",
    is.na(meta$IL4) ~ ""),
    width = unit(6.5, "mm")),
  annotation_name_gp = gpar(fontsize = 11),
  annotation_name_side = "top",
  annotation_name_rot = 0,
  show_annotation_name = TRUE
)
ha_col = HeatmapAnnotation(Col = anno_empty(which = "column", colnames(lfc_select), border = F, height = unit(5, "mm")))

p3 <- Heatmap(lfc_select,
              cell_fun = function(j, i, x, y, width, heigth, fill) {
                if(!is.na(padj_asterisk[i, j]))
                  grid.text(padj_asterisk[i,j], x, y, gp = gpar(fontsize = 10))
              },
              heatmap_legend_param = list(title = expression(paste("log" [2], "-fc"))),
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_dend = FALSE, row_names_side = "left",
              right_annotation = ha,
              top_annotation = ha_col
)

p3@right_annotation@anno_list[["LPS"]]@name_param[["show"]] <- TRUE
p3@right_annotation@anno_list[["IL4"]]@name_param[["show"]] <- TRUE

p3
