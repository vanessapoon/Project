#always setwd before starting
setwd("~/OneDrive - University of Birmingham/Year 3/Research Project/data") 

library(tidyverse)

#check files to work with
list.files("kupffer dataset/salmon") 
salmon.files <- list.files("kupffer dataset/salmon", full=TRUE)
salmon.files

read.delim("kupffer dataset/SRA.tabular") #reads files in table format
samples <- read.delim("kupffer dataset/SRA.tabular")
samples

samples %>% arrange(Run) #orders rows by Run and select columns of interest
samples <- samples %>% arrange(Run) %>% select("Run", "Tissue", "disease_state")
row.names(samples) <- samples$Run
samples

library(tximport)
library(GenomicFeatures)

#no spaces when using genomic features package so change variable type and rename
samples$Tissue <- as.factor(samples$Tissue)
samples$Tissue <- factor(samples$Tissue, labels = c("liver", "macrophage"))

gtf <- "gencode.gtf"
txdb <- makeTxDbFromGFF(gtf) #store transcript annotations
saveDb(txdb,"txdb.filename")

txdf <- AnnotationDbi::select(txdb, keys(txdb, "GENEID"), "TXNAME", "GENEID")
tab <- table(txdf$GENEID)
txdf$ntx <- tab[match(txdf$GENEID, names(tab))]
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)
txi <- tximport(salmon.files, type = "salmon", tx2gene = tx2gene, txOut = FALSE)

#groups two columns together, easier to identify sample
samples$group <- paste(samples$Tissue, samples$disease_state, sep = ".")

#differential gene expression analysis
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, samples, design = ~ group) 
dds <- DESeq(dds)
results(dds, contrast = c("group", "macrophage.Cirrhosis", "macrophage.Control")) %>%
  as.data.frame() %>%
  filter(padj<0.05)

# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")  
library(AnnotationDbi)
library(org.Hs.eg.db)

#add column with gene symbols to results data frame
res <- results(dds, contrast = c("group", "macrophage.Cirrhosis", "macrophage.Control")) %>%
  as.data.frame() %>%
  mutate(gene_id = rownames(dds)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(gene_id, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"
  )) 

# install.packages("Hmisc")
# install.packages("pheatmap")
library(Hmisc)
library(pheatmap)

#heatmap to show correlation of samples
col.data <- samples
corr.mat <- Hmisc::rcorr(counts(dds))
pheatmap::pheatmap(corr.mat$r, annotation_col = col.data)

#dispersion statistics and correction that DESeq2 applies
plotDispEsts(dds)
