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

#groups two columns together, easier to identify sample
samples$group <- paste(samples$Tissue, samples$disease_state, sep = ".")


#filtering out macrophage samples only
salmon.files <- salmon.files[samples$group%in%c("macrophage.Cirrhosis", "macrophage.Control") ]
samples <- samples[samples$group%in%c("macrophage.Cirrhosis", "macrophage.Control"),]


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

library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, samples, design = ~ group) 
dds <- DESeq(dds)
results(dds, contrast = c("group", "macrophage.Cirrhosis", "macrophage.Control")) %>%
  as.data.frame() %>%
  filter(padj<0.05)

colnames(dds) 
dds <- dds [, -c(5,7)] #excluding files 88 and 90 for downstream analysis
dds <- dds [, -c(1,4)] #excluding files 84 and 87

# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")  
library(AnnotationDbi)
library(org.Hs.eg.db)

#add column with gene symbols to results data frame
#specify name of variable, name of level for numerator and name of level for denominator, lfc= log2(cirrhosis/control)
res <- results(dds, contrast = c("group", "macrophage.Cirrhosis", "macrophage.Control")) %>%
  as.data.frame() %>%
  mutate(gene_id = rownames(dds)) %>%
  mutate(Symbol = mapIds(
    org.Hs.eg.db,
    keys = substr(gene_id, 1, 15),
    keytype = "ENSEMBL",
    column = "SYMBOL"
  )) 

#looking at multiple genes
genes <- c("GREM1", "NOG", "BMP4", "BMP7", "TGFB1", "CHRD", "BMP2")
res %>% filter(Symbol%in%genes)

#looking for one gene 
res %>% filter(Symbol=="GREM1")

# install.packages("Hmisc")
# install.packages("pheatmap")
library(Hmisc)
library(pheatmap)

#heatmap to show correlation of samples
col.data <- samples %>% dplyr::select("group")
corr.mat <- Hmisc::rcorr(counts(dds))
pheatmap::pheatmap(corr.mat$r, annotation_col = col.data)

#dispersion statistics and correction that DESeq2 applies
plotDispEsts(dds)

#counts of each sample
colSums(counts(dds))

#filtering genes from res and remove any which are NA
gene_list <- res$gene_id [res$padj < 0.05 & !is.na(res$padj)] 
gene_list <- substr(gene_list, 1, 15) 
gene_list

#BiocManager::install("clusterProfiler")
library(clusterProfiler)
enr_GO <- enrichGO(gene_list, org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
enr_GO <- setReadable(enr_GO, org.Hs.eg.db, keyType = "ENSEMBL") #converting to gene symbols
head(enr_GO, n=Inf)
