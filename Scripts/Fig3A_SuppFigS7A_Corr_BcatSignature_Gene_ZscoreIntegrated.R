library(tidyr)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
source("./Scripts/DEanalysisAndCo_Functions.R")


# 1. CORRELATION GENE-SIGNATURE AND GSEAs ------------------------------------

### Settings ----------------------------------------------------------------
Correlation_type <- "pearson"
geneSet <- "Hallmarks"
norm <- "Zscore"

signature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3", "TNFRSF19")

inpath <- "./ACC_datasets_integration/"
out <- paste0(inpath,"GSEAs/")

### Input --------------------------------------------------------------
data <- readRDS(file=paste0(inpath, "Counts_Zscore_allDatasets.RDS"))
metadata <- readRDS(file=paste0(inpath,"Metadata_allDatasets.RDS"))
# Put column of GeneSymbol as first
col_genenames <- which(colnames(data) == "GeneSymbol")
data_ok <- as.data.frame(cbind(data[,col_genenames], data[,-col_genenames]))
colnames(data_ok)[1] <- "GeneSymbol" 

### Calculate mean signature expression ---------------
avg_signature <- function(dataset, signature){
  subset <- dataset[dataset$GeneSymbol %in% signature,]
  mean <- as.data.frame(apply(subset[,-1], 2, mean))
  mean <- as.data.frame(cbind(rownames(mean), mean))
  colnames(mean) <- c("Sample","SignatureMean")
  mean_table <- as.data.frame(pivot_wider(mean, names_from=Sample, values_from=SignatureMean))
  mean_table <- cbind(c("signature_mean"), mean_table)
  colnames(mean_table)[1] <- "GeneSymbol" 
  return(mean_table)
}
signature_mean <- avg_signature(dataset=data_ok, signature=signature) # Calculate signature mean expression

### Calculate correlations between each measured gene and the signature ---------------------------------------------
all_correlations <- as.data.frame(matrix(ncol=4))
colnames(all_correlations) <- c("gene","Correlation_type","Correlation_value","pvalue")

dataset_signature <- as.data.frame(rbind(data_ok, signature_mean))
rownames(dataset_signature) <- dataset_signature[,1]
dataset_signature <- as.matrix(dataset_signature[,-1])
index_signature <- nrow(dataset_signature)
index_lastGene <- nrow(dataset_signature)-1
for(gene in 1:index_lastGene){
      tmp <- c(rownames(dataset_signature)[gene], Correlation_type, 
               round(cor.test(dataset_signature[gene,], dataset_signature[index_signature,], method=Correlation_type)$estimate, digits=2),
               cor.test(dataset_signature[gene,], dataset_signature[index_signature,], method=Correlation_type)$p.value)
      names(tmp) <- c("gene","Correlation_type","Correlation_value","pvalue")
      all_correlations <- rbind(all_correlations, tmp)
 }
all_correlations <- all_correlations[-1,]
all_correlations$BH <- p.adjust(all_correlations$pvalue, method="BH")
all_correlations$Bonferroni <- p.adjust(all_correlations$pvalue, method="bonferroni")

saveRDS(all_correlations, file=paste0(out, Correlation_type,"Correlation_BcateninActivationSignatureVsGene_WithPvalues.RDS"))


### GSEA on correlations ----------------------------------------------------
print("--------Run GSEAs---------------")
geneSet <- "Hallmarks"

all_correlations <- readRDS(file=paste0(out, Correlation_type,"Correlation_BcateninActivationSignatureVsGene_WithPvalues.RDS"))
RankedList <- as.numeric(as.vector(all_correlations$Correlation_value))
names(RankedList) <- all_correlations$gene
RankedList <- sort(RankedList, decreasing=T)  ### Rank genes by correlation
saveRDS(file=paste0(out, "RankedList_", Correlation_type,".RDS"), RankedList)
GSEA_hallmark_res <- GSEA_MSigDBHallmark_v2(bulkType=paste0("bulkIntegrated,", norm, "_", Correlation_type), 
                                            analysis=paste0("bulkIntegrated,", norm, "_", Correlation_type),
                                            RankedList=RankedList, RankedListType=paste0(Correlation_type,"Correlation"), outPath=out,
                                            geneIdType="gene_symbol")
saveRDS(file=paste0(out, "HallmarksGSEAonCorrelationToBcateninActivationSignature.RDS"), GSEA_hallmark_res)

### Plot SIGNIFICANT (padj<0.1)  pathways associated with B-catenin activation -----------
gsea <- readRDS(file=paste0(out, "/HallmarksGSEAonCorrelationToBcateninActivationSignature.RDS"))
sign <- gsea[gsea$padj < 0.1,]
sign_pos <- sign[sign$NES > 0,]
sign_pos <- as.data.frame(cbind(sign_pos, rep("pos", nrow(sign_pos)), rep(geneSet, nrow(sign_pos))))
colnames(sign_pos)[(ncol(sign_pos)-1):ncol(sign_pos)] <- c("Sign","GeneSet")
sign_neg <- sign[sign$NES < 0,]
sign_neg <- as.data.frame(cbind(sign_neg, rep("neg", nrow(sign_neg)), rep(geneSet, nrow(sign_neg))))
colnames(sign_neg)[(ncol(sign_neg)-1):ncol(sign_neg)] <- c("Sign","GeneSet")

pathways_order <- c("HALLMARK_APOPTOSIS",
                    "HALLMARK_PROTEIN_SECRETION",
                    "HALLMARK_ADIPOGENESIS",
                    "HALLMARK_COAGULATION",
                    "HALLMARK_HEME_METABOLISM",
                    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                    "HALLMARK_BILE_ACID_METABOLISM",
                    "HALLMARK_COMPLEMENT",
                    "HALLMARK_FATTY_ACID_METABOLISM",
                    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                    "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_ALLOGRAFT_REJECTION",
                    "HALLMARK_APICAL_JUNCTION",
                    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_UV_RESPONSE_DN",
                    "HALLMARK_TGF_BETA_SIGNALING", 
                    "HALLMARK_SPERMATOGENESIS",
                    "HALLMARK_ANGIOGENESIS",
                    "HALLMARK_MYC_TARGETS_V2","HALLMARK_MYC_TARGETS_V1",
                    "HALLMARK_E2F_TARGETS",
                    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                    "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT",
                    "HALLMARK_WNT_BETA_CATENIN_SIGNALING")

png(file=paste0(out, "GSEA_sign_pearson_ZscoreIntegratedData.png"))
ggplot(data=sign, aes(y=factor(pathway, levels=pathways_order), x=NES, fill=padj)) + geom_col() +
  theme_classic() + 
  ggtitle(paste0("GSEA significant (padj<0.1) enrichments\non ",Correlation_type," correlation signature-gene expression, ", norm, " integrated datsets"))
dev.off()

### Extract leading edges of pathways of interest ---------------------------
pathways <- c("HALLMARK_MYC_TARGETS_V1",
              "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_ALLOGRAFT_REJECTION",
              "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "HALLMARK_COMPLEMENT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_ANGIOGENESIS"
              )

for(p in pathways){
  genes <- unlist(gsea[gsea$pathway == p,]$leadingEdge)
  print(p)
  print(head(genes))
  write.table(file= paste0(out, "LeadingEdge_gsea_ZscoreNormIntegrated_", p, ".txt"), genes, quote= F, col.names = F, row.names = F, sep="\n")
}

### Plot GSEAs curves for pathways of interest-------------------
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)

inpathGSEAs <- paste0("./ACC_datasets_integration/GSEAs/")
pathways <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_WNT_BETA_CATENIN_SIGNALING")

for(path in pathways){
  print(path)
    
    png(paste0(inpathGSEAs, "CurvesGSEAs_",path,".png"))
    data <- readRDS(paste0(inpathGSEAs, "HallmarksGSEAonCorrelationToBcateninActivationSignature.RDS"))
    Ranks <- readRDS(file=paste0(inpathGSEAs, "RankedList_pearson.RDS"))
    
    i <- which(data$pathway==path)
    pval <- signif(data[i,]$pval, digits=3)
    padj <- signif(data[i,]$padj, digits=3)
    ES <- round(data[i,]$ES, digits=2)
    NES <- round(data[i,]$NES, digits=2)
    
    print(plotEnrichment(msigdbr_list_H_nr[[path]],
                         Ranks) + labs(title=ggtitle(paste0(path, "\n p=", pval, ", padj=", padj, ", ES=", ES, ", NES=", NES))))
    dev.off()
    
}



# 2. HEATMAP NORM EXPR GENES OF INTEREST IN INTEGARTED DATASET -----------------------------------------
# Plot expression of genes in leading edge of "Hallmarks_Inflammatory_Response" gene set (from GSEA on gene-signature correlations)
leadEdgeInflammatoryResponse <- scan(file="./ACC_datasets_integration/GSEAs/LeadingEdge_gsea_ZscoreNormIntegrated_HALLMARK_INFLAMMATORY_RESPONSE.txt", what=character())
GenesAnalyzed <- "leadEdgeInflammatoryResponse"
genes <- leadEdgeInflammatoryResponse
dataset_name <- "AllIntegratedZscore"

data <- readRDS(file="./ACC_datasets_integration/Counts_Zscore_allDatasets.RDS")
metadata <- readRDS(file="./ACC_datasets_integration//Metadata_allDatasets.RDS")
colnames(metadata)[1] <- "Sample"
metadata$BcatStatus <- metadata$BcatStatusGeneral 
metadata <- metadata[metadata$BcatStatus %in% c("Active","Inactive","ZNRF3_mut",NA),] 
metadata[metadata$BcatStatus %in% c("Active"),]$BcatStatus <- "CTNNB1mut_NuclBCat"
metadata[metadata$BcatStatus %in% c("Inactive"),]$BcatStatus <- "CTNNB1wt_MembrBCat"
metadata[metadata$BcatStatus %in% c("ZNRF3_mut"),]$BcatStatus <- "ZNRF3_altered"
metadata[is.na(metadata$BcatStatus),]$BcatStatus <- "unknown"

subset <- data[data$GeneSymbol %in% genes,]
rownames(subset) <- subset$GeneSymbol
indexGeneSymbol <- which(colnames(data) == "GeneSymbol")
subset <- subset[,-indexGeneSymbol]
subset <- as.matrix(subset)

# Build annotation object
all_anno <- metadata[,c('Sample','BcatStatus','BcatSignatureMeanZScore')]
all_anno <- unique(all_anno) # required to remove redundancy in tcga where metadata for TCGA-OR-A5KB has 2 lines beacuse 2 different mutations on TP53
rownames(all_anno) <- all_anno$Sample
all_anno$BcatStatus <- factor(all_anno$BcatStatus, levels=c("CTNNB1mut_NuclBCat","CTNNB1wt_MembrBCat","ZNRF3_altered", "unknown"))
samples <- intersect(colnames(subset), all_anno$Sample) ### select only samples present both in matrix of counts and annotation: note that samples with BcatStatus unknown or ZNRF3mut are remooved
subset_order <- subset[,match(samples, colnames(subset))] ### set same order for matrix of counts and annotation
all_anno <- all_anno[match(samples, all_anno$Sample),]
subset_order <- as.matrix(subset_order) ### Transform input to heatmap into matrix, as required

col_fun = colorRamp2(c(-1.7, 0,1.7), c("green", "white", "red"))
annot <- HeatmapAnnotation(BcatStatus=all_anno$BcatStatus,
                           BcatSignatureMeanZScore=all_anno$BcatSignatureMeanZScore,
                           col=list(BcatStatus=c("CTNNB1mut_NuclBCat"="red","ZNRF3_altered"="orange","unknown"="grey","CTNNB1wt_MembrBCat"="green"),
                                    BcatSignatureMeanZScore=col_fun))
png(paste0(out, "HeatmapGeneExpression_", GenesAnalyzed, dataset_name, ".png"), height = 1200, width=800)
print(ComplexHeatmap::Heatmap(subset_order
                              ,column_km=2
                              ,column_title = paste0(dataset_name, "_", GenesAnalyzed) # title
                              ,top_annotation = annot)) # set sample annotation
dev.off()


# 3. GSEA IMMUNE RELATED GENE SETS ON CTNNB1 MUT vs WT FOLD CHANGES  ------------------------
Correlation_type <- "pearson"
out <- "./DEanalyses_ACC/CTNNB1mutVsWt/GSEAs/"

#Assie <- readRDS("./ACC_datasets_formatted/GSEAs_onDEanalyses/AssieFoldChange_BcatActive_vs_not_Hallmarks.RDS")  
Assie <- readRDS("./DEanalyses_ACC/CTNNB1mutVsWt/GSEAs/AssieFoldChange_BcatActive_vs_not_Hallmarks.RDS")  
Assie$Dataset <- "Assie"
#Heaton <- readRDS("./ACC_datasets_formatted/GSEAs_onDEanalyses//HeatonFoldChange_BcatActive_vs_not_Hallmarks.RDS")  
Heaton <- readRDS("./DEanalyses_ACC/CTNNB1mutVsWt/GSEAs/HeatonFoldChange_BcatActive_vs_not_Hallmarks.RDS")  
Heaton$Dataset <- "Heaton"
#tcga <- readRDS("./ACC_datasets_formatted/GSEAs_onDEanalyses/TCGAFoldChange_BcatActive_vs_not_Hallmarks.RDS")
tcga <- readRDS("./DEanalyses_ACC/CTNNB1mutVsWt/GSEAs/tcgaFoldChange_BcatActive_vs_not_Hallmarks.RDS")  
tcga$Dataset <- "tcga"

input <- rbind(Assie, Heaton, tcga)
pathways <- c("HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_IL6_JAK_STAT3_SIGNALING","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_IL2_STAT5_SIGNALING")
input <- input[input$pathway %in% pathways,]

all <- input[,c("pathway","pval","padj","log2err","ES","NES","size","leadingEdge","Dataset")]
all$pathway <- factor(all$pathway, levels=pathways)

all$significant <- ifelse(all$padj < 0.1, "*", "")
png(paste0(out, "ImmuneRelatedPathways_Assie_Heaton_tcga.png"), width=1400, height = 150)
ggplot(data=all, aes(x=NES, y=factor(Dataset, levels=c("tcga","Heaton","Assie")), fill=padj)) +
  geom_col() + ylab("dataset") + theme_bw() + facet_wrap(.~pathway, ncol = 6) +
  ggtitle("Significant (FDR<0.1) immune-related pathways")
dev.off()
