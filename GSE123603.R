setwd("~/OneDrive - University of Birmingham/Year 3/Research Project/data") 

library(tidyverse)

list.files ("GSE123603_IL4IL6/salmon") 
salmon.files <- list.files ("GSE123603_IL4IL6/salmon", full=TRUE)
salmon.files

read.delim ("GSE123603_IL4IL6/SRA.tabular")
samples <- read.delim("GSE123603_IL4IL6/SRA.tabular")
samples

samples %>% arrange(Run) 
samples <- samples %>% arrange(Run) %>% dplyr::select("Run", "Cell_type", "Treatment") 
row.names(samples) <- samples$Run
samples

samples$Treatment <- c("control", "control", "control", "IL4", "IL4", "IL4", "IL6", "IL6", "IL6", "both", "both", "both")

samples$Replicate <- c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3") #to allow for batch effect correction
samples

library(tximport)
library(GenomicFeatures)

samples$Cell_type <- as.factor(samples$Cell_type)
samples$Cell_type <- factor(samples$Cell_type, labels = c("macrophage"))

samples$group <- paste(samples$Cell_type, samples$Treatment, sep = ".")
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
dds$group <- relevel(dds$group, ref ="macrophage.control") 
dds <- DESeq(dds)
results(dds, contrast = c("group", "macrophage.IL4", "macrophage.control")) %>%  #change to look at IL4, IL6 or both
  as.data.frame() %>%
  filter(padj<0.05)

library(AnnotationDbi)
library(org.Hs.eg.db)

res <- results(dds, contrast = c("group", "macrophage.IL4", "macrophage.control")) %>%  #change to look at IL4, IL6 or both
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

genes <- c("CD14", "ITGAM", "FCGR3A", "CD163", "CD3D", "CD19", "PTPRC")
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

dds$group <- relevel(dds$group, ref = "macrophage.control")

resultsNames(dds)

library (apeglm)

res <- lfcShrink(dds, coef =  "group_macrophage.IL4_vs_macrophage.control") %>% as.data.frame()

results_all <- readRDS("results_all.RDS")
results_all[["GSE123603_IL4"]] <- res

saveRDS(results_all, file = "results_all.RDS")

results_all


