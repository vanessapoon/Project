setwd("~/OneDrive - University of Birmingham/Year 3/Research Project/data") 

library(tidyverse)

list.files ("GSE100382_LPSTNFIFNa/Salmon") 
salmon.files <- list.files ("GSE100382_LPSTNFIFNa/Salmon", full=TRUE)
salmon.files

read.csv ("GSE100382_LPSTNFIFNa/samples.csv")
samples <- read.csv ("GSE100382_LPSTNFIFNa/samples.csv")
samples

samples %>% arrange(Sample.Name) 
row.names(samples) <- samples$Sample.Name
samples

library(tximport)
library(GenomicFeatures)

samples$group <- paste(samples$Cell_Type, samples$Treatment, sep = ".")
samples

gtf <- "gencode.gtf"
txdb <- makeTxDbFromGFF(gtf)
saveDb(txdb,"txdb.filename")


txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)
txi <- tximport(salmon.files, type = "salmon", tx2gene = tx2gene, txOut = FALSE)

library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, samples, design = ~ group) 
dds$group <- relevel(dds$group, ref ="macrophage.None") 
dds <- DESeq(dds)
results(dds, contrast = c("group", "macrophage.IFNLPSTNF", "macrophage.None")) %>%  #change to look at separate treatments
  as.data.frame() %>%
  filter(padj<0.05)

library(AnnotationDbi)
library(org.Hs.eg.db)

res <- results(dds, contrast = c("group", "macrophage.IFNLPSTNF", "macrophage.None")) %>%  #change to look at separate treatments
  as.data.frame() %>%
  mutate(gene_id = rownames(dds)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(gene_id, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"
  )) 


genes <- c("IL6", "IL12A", "TNF", "IL10", "ARG1", "CCL17", "CCL22")
res %>% filter(Symbol%in%genes)

res %>% filter(Symbol=="GREM1")

library(Hmisc)
library(pheatmap)

col.data <- samples %>% dplyr::select("group")
corr.mat <- Hmisc::rcorr(counts(dds))
pheatmap::pheatmap(corr.mat$r, annotation_col = col.data)

plotDispEsts(dds)
colSums(counts(dds))

gene_list <- res$gene_id [res$padj < 0.05 & !is.na(res$padj)] 
gene_list <- substr(gene_list, 1, 15) 
gene_list

library(clusterProfiler)
enr_GO <- enrichGO(gene_list, org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
enr_GO <- setReadable(enr_GO, org.Hs.eg.db, keyType = "ENSEMBL") #converting to gene symbols
head(enr_GO, n=Inf)

resultsNames(dds)

res <- lfcShrink(dds, coef =  "group_macrophage.IFNLPSTNF_vs_macrophage.None") %>% as.data.frame()

#results_all <- readRDS("results_all.RDS")
#results_all[["GSE100382_LPS"]] <- res

#saveRDS(results_all, file = "results_all.RDS")


#saving other results in separate file 

#results_GSE100382 <- list()
#saveRDS(results_GSE100382, file = "results_GSE100382.RDS")

results_GSE100382 <- readRDS("results_GSE100382.RDS")
results_GSE100382[["GSE100382_IFNLPSTNF"]] <- res

saveRDS(results_GSE100382, file = "results_GSE100382.RDS")


library(cowplot)
library(hrbrthemes)

results_GSE100382 <- readRDS("results_GSE100382.RDS")

names(results_GSE100382) <- plyr::mapvalues(x = names(results_GSE100382), from = "GSE100382_IFNLPSTNF", to = "GSE100382_IFNTNFLPS")

data <- read.delim("GSE100382.txt")

sub <- data$ID

lfc <- list()
for (i in sub){
  lfc[[i]] <- results_GSE100382[[i]]$log2FoldChange
}

padj <- list()
for (i in sub){
  padj[[i]] <- results_GSE100382[[i]]$padj
}

genes <- rownames(results_GSE100382[[1]])

lfc_df <- lfc %>% as.data.frame()
rownames(lfc_df) <- genes

padj_df <- padj %>% as.data.frame()
rownames(padj_df) <- genes

library(org.Hs.eg.db)
Symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = substr(genes, 1, 15), keytype = "ENSEMBL", column = "SYMBOL")

targets <- c("GREM1")

lfc_select <- lfc_df [Symbol %in% targets, ]
rownames(lfc_select) <- Symbol [Symbol %in% targets]
lfc_long <- lfc_select %>% mutate(gene_id = factor(rownames(lfc_select))) %>% pivot_longer(cols = 1:7)
padj_select <- padj_df [Symbol %in% targets, ]
rownames(padj_select) <- Symbol [Symbol %in% targets]
padj_long <- padj_select %>% mutate(gene_id = rownames(padj_select)) %>% pivot_longer(cols = 1:7)
padj_long$value [is.na(padj_long$value)] <- 1

lfc_long$padj <- padj_long$value

lfc_long$signif <- case_when(
  lfc_long$padj < 0.001 ~ "***",
  lfc_long$padj < 0.01 ~ "**",
  lfc_long$padj < 0.05 ~ "*") 

max <- max(abs(lfc_long$value), na.rm = T)
p2 <- lfc_long %>% 
  ggplot(aes(x = factor(name, level = c("GSE100382_LPS", "GSE100382_TNF", "GSE100382_IFN", "GSE100382_TNFLPS", "GSE100382_IFNLPS", "GSE100382_IFNTNF", "GSE100382_IFNTNFLPS")), y = gene_id, fill = value))+
  geom_tile()+
  geom_text(aes(label = signif))+
  labs(title = "", x ="", y = "", 
       fill = expression(paste("log" [2], "-fc")))+
  scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-max, max)) +
  scale_y_discrete(limits = rev(levels(lfc_long$gene_id)))+
  theme_minimal()+
  theme(line = element_line(colour = "black", size = 0.5),
        # text = element_text(family = "serif"),
        # plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.background = element_blank(), panel.grid = element_blank(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.key = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.ticks = element_line(size = 0.5, colour = "black"),
        axis.title = element_text(size = 10, hjust = 0.5),
        axis.text.x = element_text(size=10, angle = 90, hjust = 1),
        axis.text.y = element_text(size=10, ))
p2

#compelx heatmap instead
lfc_select <- lfc_df [Symbol %in% targets, ] %>% as.matrix()
rownames(lfc_select) <- Symbol [Symbol %in% targets]
padj_select <- padj_df [Symbol %in% targets, ] %>% as.matrix()
rownames(padj_select) <- Symbol [Symbol %in% targets]


padj_asterisk <- matrix(case_when(
  padj_select < 0.001 ~ "***",
  padj_select < 0.01 ~ "**",
  padj_select < 0.05 ~ "*"), ncol = ncol(padj_select))

p4 <- Heatmap(lfc_select,
              cell_fun = function(j, i, x, y, width, heigth, fill) {
                if(!is.na(padj_asterisk[i, j]))
                  grid.text(padj_asterisk[i,j], x, y, gp = gpar(fontsize = 10))
              },
              heatmap_legend_param = list(title = expression(paste("log" [2], "-fc"))),
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_dend = FALSE, row_names_side = "left",
              
)
p4

