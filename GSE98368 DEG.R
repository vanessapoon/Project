setwd("~/OneDrive - University of Birmingham/Year 3/Research Project/data") 

library(tidyverse)

list.files ("GSE98368_IFNg/Salmon") 
salmon.files <- list.files ("GSE98368_IFNg/Salmon", full=TRUE)
salmon.files

read.delim ("GSE98368_IFNg/SRA.tabular")
samples <- read.delim("GSE98368_IFNg/SRA.tabular")
samples

samples %>% arrange(Run) 
samples <- samples %>% arrange(Run) %>% select("Run", "Cell_type", "Treatment") 
row.names(samples) <- samples$Run
samples

library(tximport)
library(GenomicFeatures)

samples$Cell_type <- as.factor(samples$Cell_type)
samples$Cell_type <- factor(samples$Cell_type, labels = c("macrophage"))

samples$group <- paste(samples$Cell_type, samples$Treatment, sep = ".")

salmon.files <- salmon.files[samples$group%in%c("macrophage.IFN-?", "macrophage.none") ]
samples <- samples[samples$group%in%c("macrophage.IFN-?", "macrophage.none"),]

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
dds <- DESeq(dds)
results(dds, contrast = c("group", "macrophage.IFN-?", "macrophage.none")) %>%
  as.data.frame() %>%
  filter(padj<0.05)

colnames(dds) 

library(AnnotationDbi)
library(org.Hs.eg.db)

res <- results(dds, contrast = c("group", "macrophage.IFN-?", "macrophage.none")) %>%
  as.data.frame() %>%
  mutate(gene_id = rownames(dds)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(gene_id, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"
  )) 

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

