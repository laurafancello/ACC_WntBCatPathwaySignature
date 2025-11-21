library("SummarizedExperiment")
library("DESeq2")
library("IHW") 
library("biomaRt")
library("apeglm") 
library("pheatmap") 
library("RColorBrewer")
library("reshape2")
library("stringr")
library("TCGAretriever")
library("ggplot2")
library("ggrepel")
library("dplyr")

outPath <- "./DEanalyses_ACC/CTNNB1mutVsWt/"

# Read input raw counts ---------------------------------------------------
rna <- read.table(file="./ACC_datasets_formatted/tcga/RNAseq_STAR_RawCounts_fromTCGABiolinks_GeneSymbols.txt", sep="\t", header=T)
colnames(rna) <- str_replace_all(colnames(rna), '\\.', "-")

CTNNB1mutWt <- readRDS("./ACC_datasets_formatted/tcga/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS") # Read manually curated file with CTNNB1 or ZNRF3 altered samples 
rna <- rna[,which(colnames(rna) %in% CTNNB1mutWt$Sample)] # keep only samples with both RNA and MUT/CN information
CTNNB1mutWt <- CTNNB1mutWt[CTNNB1mutWt$Alteration %in% c("CTNNB1_mut","wt"),] # keep only CTNNB1mut and wt samples, discard ZNRF3 altered samples
CTNNB1mutWt$Alteration <- factor(CTNNB1mutWt$Alteration, levels=c("wt","CTNNB1_mut"))

# Order samples in counts and metadata the same way
rna <- rna[,match(CTNNB1mutWt$Sample,colnames(rna))]

if(identical(CTNNB1mutWt$Sample, colnames(rna))){
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = rna,
    colData = CTNNB1mutWt,
    design = ~ Alteration)
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds) 
  res <- results(dds)
  res <- res[!(is.na(res$padj)),]
  write.csv(file=paste0(outPath,"DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_tcga.csv"), res)
}
