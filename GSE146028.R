setwd("~/OneDrive - University of Birmingham/Year 3/Research Project/data") 

library(tidyverse)

list.files ("GSE146028_IFNgLPSIL4/salmon") 
salmon.files <- list.files ("GSE146028_IFNgLPSIL4/salmon", full=TRUE)
salmon.files

read.delim ("GSE146028_IFNgLPSIL4/SRA.tabular")
samples <- read.delim("GSE146028_IFNgLPSIL4/SRA.tabular")
samples

samples %>% arrange(Run) 
samples <- samples %>% arrange(Run) %>% dplyr::select("Run", "Cell_type", "treatment.cytokine") 
row.names(samples) <- samples$Run
samples

library(tximport)
library(GenomicFeatures)

samples$Cell_type <- as.factor(samples$Cell_type)
samples$Cell_type <- factor(samples$Cell_type, labels = c("macrophage"))

samples$group <- paste(samples$Cell_type, samples$treatment.cytokine, sep = ".")
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
dds$group <- relevel(dds$group, ref ="macrophage.none") 
dds <- DESeq(dds)
results(dds, contrast = c("group", "macrophage.IL4", "macrophage.none")) %>%  #change to look at separate treatments
  as.data.frame() %>%
  filter(padj<0.05)

library(AnnotationDbi)
library(org.Hs.eg.db)

res <- results(dds, contrast = c("group", "macrophage.IL4", "macrophage.none")) %>%  #change to look at separate treatments
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

res <- lfcShrink(dds, coef =  "group_macrophage.IL4_vs_macrophage.none") %>% as.data.frame()

results_all <- readRDS("results_all.RDS")
results_all[["GSE146028_IL4"]] <- res

saveRDS(results_all, file = "results_all.RDS")

names(results_all) <- plyr::mapvalues(x = names(results_all), from = "GSE146028_IFNg", to = "GSE146028_IFNg+LPS")

