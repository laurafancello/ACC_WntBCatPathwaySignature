library(tidyr)
library(ggplot2)
library(stats)
library(stringr)
library(gridExtra)
library(ComplexHeatmap)

out <- "./SignatureAndLEF1_inACC/output/"
analysis <- "BcatActive_Inactive"

# Read genes of interest: known B-catenin targets ------------------------
knownTargets <- c("ABCB1", "AFF3", "AXIN2", "BCL2L2", "BIRC5", "CCND1", "CDC25A", "CDKN2A", "CDX1", "CLDN1", "CTLA4", "DKK1", "EDN1", "ENAH", "ENC1", "FGF18", 
                 "FGF4", "FGFBP1", "FOSL1", "FSCN1", "FST", "FZD7", "GBX2", "HES1", "HNF1A", "ID2", "JAG1", "JUN", "KRT5", "L1CAM","LAMC2", "LEF1", "LGR5", "MMP14", 
                 "MMP7", "MYC", "MYCBP", "NEDD9", "NEUROD1", "NEUROG1", "NOS2", "NOTCH2", "NRCAM", "PDE2A", "PLAU", 
                 "PPARD", "PTGS2", "S100A4", "SGK1", "SMC3", "SP5", "SUZ12", "TBX3", "TBXT", "TCF4", "TERT", "TNC", "TNFRSF19", "VCAN", "VEGFA", 
                 "YY1AP1") # known B-catenin targets (those reported on Nusse lab website + Herbst et al 2014 paper + PDE2A + AFF3)
genes_of_interest <- list(knownTargets=knownTargets)

# Read input log2 norm counts  --------
### For microarray data, if multiple probes for one gene use their average
in_Assie <- readRDS(file="./ACC_datasets_formatted/Assie/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")
in_Heaton <- readRDS(file="./ACC_datasets_formatted/Heaton/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")
in_tcga <- readRDS(file="./ACC_datasets_formatted/tcga/Log2Pseudocount1_DESeqNormCounts_tcga.RDS")
in_Lefevre <- readRDS(file="./ACC_datasets_formatted/Lefevre/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")

# Read sample classification and filter----------------------------------------------
samples_BcatMutWt_TCGA <- readRDS(file="./ACC_datasets_formatted/tcga/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
samples_BcatMutWt_TCGA <- samples_BcatMutWt_TCGA[samples_BcatMutWt_TCGA$Alteration %in% c("CTNNB1_mut","wt"),]
samples_BcatMutWt_Assie <- as.data.frame(readRDS(file="./ACC_datasets_formatted/Assie/Metadata/Assie_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS"))
samples_BcatMutWt_Assie <- samples_BcatMutWt_Assie[samples_BcatMutWt_Assie$Alteration %in% c("CTNNB1","wt"),]
samples_BcatMutWt_Heaton <- readRDS(file="./ACC_datasets_formatted/Heaton/Metadata/Heaton_sampleIDs_BcatNuclearMembrane_onlyACCs.RDS")
samples_BcatMutWt_Heaton <- samples_BcatMutWt_Heaton[samples_BcatMutWt_Heaton$Alteration %in% c("Membrane","Nuclear"),]
meta_Lefevre <- read.csv(file="./ACC_datasets_formatted/Lefevre/Metadata/Metadata_Lefevre.csv")
samples_BcatMutWt_Lefevre <- meta_Lefevre[!(meta_Lefevre$BcatStatus %in% c("ctrl","ctrlDoxy")),c('Assay.Name','BcatStatus')]
colnames(samples_BcatMutWt_Lefevre) <- c("Sample","Alteration")
  
### Discard ZNRF3 altered samples, discard non-ACC histotype samples
in_tcga_select <- in_tcga[,c(1,which(colnames(in_tcga) %in% samples_BcatMutWt_TCGA$Sample))]
in_Assie_select <- in_Assie[,c(1,which(colnames(in_Assie) %in% samples_BcatMutWt_Assie$Sample))]
in_Heaton_select <- in_Heaton[,c(1,which(colnames(in_Heaton) %in% samples_BcatMutWt_Heaton$Sample))]

Assie_levels <- c("CTNNB1","wt")
Heaton_levels <- c("Nuclear","Membrane")
tcga_levels <- c("CTNNB1_mut","wt")

# Generate input heatmaps
in_Assie_long  <- pivot_longer(in_Assie_select, cols=2:ncol(in_Assie_select), names_to="Sample", values_to="Value")
input_boxplot_Assie <- merge(in_Assie_long,samples_BcatMutWt_Assie,by="Sample",all=F)
input_boxplot_Assie <- input_boxplot_Assie[!(is.na(input_boxplot_Assie$Alteration)),]
input_boxplot_Assie$Value <- as.numeric(as.vector(input_boxplot_Assie$Value))

in_tcga_long <- pivot_longer(in_tcga_select, cols=2:ncol(in_tcga_select), names_to="Sample", values_to="Value")
input_boxplot_tcga <- merge(in_tcga_long,samples_BcatMutWt_TCGA,by="Sample",all=F)
input_boxplot_tcga <- input_boxplot_tcga[!(is.na(input_boxplot_tcga$Alteration)),]
input_boxplot_tcga$Value <- as.numeric(as.vector(input_boxplot_tcga$Value))

in_Heaton_long  <- pivot_longer(in_Heaton_select, cols=2:ncol(in_Heaton_select), names_to="Sample", values_to="Value")
input_boxplot_Heaton <- merge(in_Heaton_long,samples_BcatMutWt_Heaton,by="Sample",all=F)
input_boxplot_Heaton <- input_boxplot_Heaton[!(is.na(input_boxplot_Heaton$Alteration)),]
input_boxplot_Heaton$Value <- as.numeric(as.vector(input_boxplot_Heaton$Value))

in_Lefevre_long  <- pivot_longer(in_Lefevre, cols=2:ncol(in_Lefevre), names_to="Sample", values_to="Value")
input_boxplot_Lefevre <- merge(in_Lefevre_long, samples_BcatMutWt_Lefevre, by="Sample",all=F)
input_boxplot_Lefevre$Value <- as.numeric(as.vector(input_boxplot_Lefevre$Value))

for(i in 1:length(genes_of_interest)){
  GenesAnalyzed <- names(genes_of_interest)[i]
  genes <- genes_of_interest[[i]]
  print(GenesAnalyzed)

  # Collect Wilcoxon p-values
  all <- as.data.frame(matrix(ncol=7))
    
  for(dataset in c("Assie","Heaton","tcga")){
    
    print(dataset)
    
    if(dataset == "Assie"){
      input_boxplot <- input_boxplot_Assie
      levels <- Assie_levels
    }
    if(dataset == "Heaton"){
      input_boxplot <- input_boxplot_Heaton
      levels <- Heaton_levels
    }
    if(dataset == "tcga"){
      input_boxplot <- input_boxplot_tcga
      levels <- tcga_levels
    }
    
    input_boxplot$Value <- as.numeric(as.vector(input_boxplot$Value))
    input_boxplot$Alteration <- factor(input_boxplot$Alteration, levels=levels)
    
    for(gene in genes){
      print(gene)
      
      if(nrow(input_boxplot[input_boxplot$GeneSymbol==gene,])==0){
        all <- rbind(all, c(gene, GenesAnalyzed, analysis, dataset, "rawP_Wilcoxon", NA, NA))
      }else{
        df_sub <- input_boxplot[((input_boxplot$GeneSymbol==gene)&(input_boxplot$Alteration %in% levels[1:2])),]
        pval <- signif(wilcox.test(Value ~ Alteration, data = df_sub)$p.value, 3)
        FC <- mean(df_sub[df_sub$Alteration %in% levels[1],]$Value) - mean(df_sub[df_sub$Alteration %in% levels[2],]$Value)
        all <- rbind(all, c(gene, GenesAnalyzed, analysis, dataset, "rawP_Wilcoxon", pval, FC))
      }
    }
  }
}

dataset <- "Lefevre"
print(dataset)
levels <- c("expressed", "repressed")

for(gene in genes){
  print(gene)
  if(nrow(input_boxplot_Lefevre[input_boxplot_Lefevre$GeneSymbol==gene,])==0){
    all <- rbind(all, c(gene, GenesAnalyzed, analysis, dataset, "rawP_Ttest", NA, NA))
  }else{
    df_sub <- input_boxplot_Lefevre[((input_boxplot_Lefevre$GeneSymbol==gene)&(input_boxplot_Lefevre$Alteration %in% levels[1:2])),]
    df_sub$Value <- as.numeric(as.vector(df_sub$Value))
    pval <- signif(t.test(Value ~ Alteration, data = df_sub)$p.value, 3)
    FC <- mean(df_sub[df_sub$Alteration %in% levels[1],]$Value) - mean(df_sub[df_sub$Alteration %in% levels[2],]$Value)
    all <- rbind(all, c(gene, GenesAnalyzed, analysis, dataset, "rawP_Ttest", pval, FC))
  }
}

all <- all[-1,]
colnames(all) <- c("Gene", "GenesAnalyzed", "analysis_type", "Dataset","Test", "Sign_value", "log2FC")
all$Sign_value <- as.numeric(as.vector(all$Sign_value))
all$log2FC <- as.numeric(as.vector(all$log2FC)) 
#saveRDS(file=paste0(out, "WilcoxonTtests_Pval_FC_",GenesAnalyzed,"_", analysis, ".RDS"), all)

##################################################################################
# Heatmap of all FC, highlighting by "*" raw P <0.05 FOR FIGURE 1B --------------------------------------------
##################################################################################
all$significant <- ifelse(all$Sign_value<0.05, "*", "")
data <- all[!(is.na(all$Gene)),]

png(file=paste0(out, "HeatmapAllFC_HighlightSigntPWilcoxOrTtest_", GenesAnalyzed, "_allPublicACCDatasets_",analysis,".png"), height = 900)
print(ggplot(data=data, aes(x=Dataset, y=Gene)) + 
        geom_tile(aes(fill=log2FC), colour="black") + 
        geom_text(aes(label=significant)) +
        ggtitle(paste0("Significant Wilcoxon (or T) test (<0.05), ", GenesAnalyzed), subtitle = analysis) +
        scale_fill_gradient2(low = "green", mid="white", high = "red") +
        theme_bw())
dev.off()

