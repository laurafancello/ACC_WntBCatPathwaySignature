library(fgsea)
library(ggplot2)
library(writexl)
library(msigdbr)
library(fgsea)
library(gridExtra)
library(ggplot2)
options(digits=15)
source("./Scripts/Miscellaneous_Functions.R")

analysis = "CTNNB1_vs_wt"

inpath <- "./DEanalyses_ACC/CTNNB1mutVsWt/"
out <- "./DEanalyses_ACC/CTNNB1mutVsWt/GSEAs/"
if (!(file.exists(out))){
  dir.create(file.path(out))
}

# 1. GENERATE RANKED LISTS FOR GSEAs  ------------------------
colNameGenes <- "ID"

for(dataset_name in c("Assie","Lefevre","Heaton", "tcga")){
  print(dataset_name)
  
  if(dataset_name == "Assie"){
    res <- read.csv(file=paste0(inpath, "./DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_Assie.csv"))
    res$log2FoldChange <- res$logFC 
    
    RankedList <- res$log2FoldChange
    names(RankedList) <- res[,colNameGenes]
    if(length( which(is.na(RankedList)))>0){
      RankedList <- RankedList[-(which(is.na(RankedList)))]
    }
    RankedList <- sort(RankedList)
    RankedList <- as.data.frame(cbind(names(RankedList),RankedList))
    colnames(RankedList) <- c("GeneName","RankedList")
    write.table(file=paste0(out, "FoldChangeRanked_", analysis, "_",dataset_name,".rnk"), RankedList, col.names = T, row.names = T, sep="\t")
    
  }

  if(dataset_name == "Lefevre"){
    # Note that for Lefevre dataset the comparison is B-Catenin inactivated by shRNA versus constitutively expressed
    res <- read.table(file="./ACC_datasets_formatted/Lefevre/Metadata/SuppTable1_DEanalysis_FC_adjPvalue.txt", header=T, sep="\t")
    colnames(res) <- c("ID", "log2FoldChange", "adjP")

    RankedList <- res$log2FoldChange
    names(RankedList) <- res[,colNameGenes]
    if(length( which(is.na(RankedList)))>0){
      RankedList <- RankedList[-(which(is.na(RankedList)))]
    }
    RankedList <- sort(RankedList)
    RankedList <- as.data.frame(cbind(names(RankedList),RankedList))
    colnames(RankedList) <- c("GeneName","RankedList")
    write.table(file=paste0(out, "FoldChangeRanked_", analysis, "_",dataset_name,".rnk"), RankedList, col.names = T, row.names = T, sep="\t")
    
  }
 
  if(dataset_name == "Heaton"){
    res <- read.csv(file=paste0(inpath, "DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_Heaton_ACCs.csv"))
    res$log2FoldChange <- res$logFC 
    
    RankedList <- res$logFC
    names(RankedList) <- res[,colNameGenes]
    if(length(which(is.na(RankedList)))>0){
      RankedList <- RankedList[-(which(is.na(RankedList)))]
    }
    RankedList <- sort(RankedList)
    RankedList <- as.data.frame(cbind(names(RankedList),RankedList))
    colnames(RankedList) <- c("GeneName","RankedList")
    write.table(file=paste0(out, "FoldChangeRanked_", analysis, "_",dataset_name,".rnk"), RankedList, col.names = T, row.names = T, sep="\t")
    
  }
  
  if(dataset_name == "tcga"){
    res <- read.csv(file=paste0(inpath, "DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_tcga.csv"))
    colnames(res) <- c("ID", "baseMean","log2FoldChange","lfcSE","stat","p.value", "padj")

    RankedList <- res$log2FoldChange
    names(RankedList) <- res[,colNameGenes]
    if(length( which(is.na(RankedList)))>0){
      RankedList <- RankedList[-(which(is.na(RankedList)))]
    }
    RankedList <- sort(RankedList)
    RankedList <- as.data.frame(cbind(names(RankedList),RankedList))
    colnames(RankedList) <- c("GeneName","RankedList")
    write.table(file=paste0(out, "FoldChangeRanked_", analysis, "_",dataset_name,".rnk"), RankedList, col.names = T, row.names = T, sep="\t")
    
  }
}


# 2. PERFORM GSEAs  ------------------------
# Enrichment analysis based on log2FC of DE annalysis CTNNB1_mut vs wt (Assie and TCGA), Bcat nuclear vs membrane (Heaton), Bcat expressed or repressed in H295R cel line (Lefevre)
RankedListType  <- "FoldChange"
datasets <- c("Assie","Lefevre","Heaton","tcga")
for(dataset_name in datasets){
    print(dataset_name)

    RankedList <- read.table(file=paste0(out, RankedListType, "Ranked_", analysis, "_",dataset_name,".rnk"), sep="\t", header=T)
    values <- RankedList$RankedList
    names(values) <- RankedList$GeneName
    RankedList <- values
    # Transfrom -Inf into zeros (when calculating log2FC no pseudocount was used)
    RankedList[which(RankedList== -Inf)] <- 0
    
    print("Running GSEA on Hallmarks")
    outHallmarks <- GSEA_MSigDBHallmark_v2(bulkType=dataset_name, analysis="CTNNB1_vs_wt", RankedList=RankedList, RankedListType=RankedListType, outPath=out, geneIdType="gene_symbol")
    saveRDS(file=paste0(out, dataset_name, RankedListType, "_BcatActive_vs_not_Hallmarks.RDS"), outHallmarks[order(outHallmarks$padj),])

}


# 3. PLOT MYC TARGETS GSEA ENRICHMENTS --------------------------------
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)

datasets <- c("Assie","Heaton","Lefevre","tcga")
pathways <- c("HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2")
RankedListType <- "FoldChange"
for(dataset_name in datasets){
  res  <- readRDS(file=paste0(out, dataset_name, "FoldChange_BcatActive_vs_not_Hallmarks.RDS"))
  Ranks <- read.table(file=paste0(out, RankedListType, "Ranked_CTNNB1_vs_wt_", dataset_name, ".rnk"), header = T, sep="\t")
  Ranks_ok <- as.numeric(as.vector(Ranks$RankedList))
  names(Ranks_ok) <- as.character(as.vector(Ranks$GeneName))
  for(path in pathways){
    sub <- res[res$pathway==path,]
    pval <- signif(sub$pval, digits=3)
    padj <- signif(sub$padj, digits=3)
    ES <- round(sub$ES, digits=2)
    NES <- round(sub$NES, digits=2)
    print(sub)
    png(file=paste0(out, "PlotGSEA_",path,"_",dataset_name,".png"))
    print(plotEnrichment(msigdbr_list_H_nr[[path]],
                       Ranks_ok) + labs(title=ggtitle(paste0(path, ", ", dataset_name, "\n p=", pval, ", padj=", padj, ", ES=", ES, ", NES=", NES))))
    dev.off()
  }
}


# 4. PLOT TOGETHER GSEAs IMMUNE RELATED GENE SETS ON CTNNB1 MUT vs WT FOLD CHANGES  ------------------------
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
