library(ggplot2)
library(ggpubr)
library(ggExtra)
library(gridExtra)
library(cutpointr)
library(ggrepel)
library(survival)
library(survminer)
library(ggsurvfit)

out <- "./SignatureAndLEF1_inACC/output/"

source("./Scripts/Miscellaneous_Functions.R")

#  Read signatures ---------------------------------------------
signature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19") # home-made signature of Wnt/B-Catenin pathway activation in ACC
Nusse_PDE2A <- c("ABCB1", "AFF3", "AXIN2", "BCL2L2", "BIRC5", "CCND1", "CDC25A", "CDKN2A", "CDX1", "CLDN1", "CTLA4", "DKK1", "EDN1", "ENAH", "ENC1", "FGF18", 
                 "FGF4", "FGFBP1", "FOSL1", "FSCN1", "FST", "FZD7", "GBX2", "HES1", "HNF1A", "ID2", "JAG1", "JUN", "KRT5", "L1CAM","LAMC2", "LEF1", "LGR5", "MMP14", 
                 "MMP7", "MYC", "MYCBP", "NEDD9", "NEUROD1", "NEUROG1", "NOS2", "NOTCH2", "NRCAM", "PDE2A", "PLAU", 
                 "PPARD", "PTGS2", "S100A4", "SGK1", "SMC3", "SP5", "SUZ12", "TBX3", "TBXT", "TCF4", "TERT", "TNC", "TNFRSF19", "VCAN", "VEGFA", 
                 "YY1AP1") # known B-catenin targets (those reported on Nusse lab website + Herbst et al 2014 paper + PDE2A + AFF3)
KeggWntPathwayGenes <- read.table(file="./SignatureAndLEF1_inACC/input/MSigDB_C2CPCurated_KEGG_WNT_SIGNALING_PATHWAY.v2024.1.Hs.txt", header=T, sep="\t") # KEGG genes for Wnt signaling pathway
KeggWntPathwayGenes <- KeggWntPathwayGenes$Gene_symbol
pancancer3gene <- c("AXIN2","SP5","TNFRSF19") # home-made 3-gene signature of Wnt/B-Catenin pathway activation for multiple cancer types
                 
genesOfInterest <- c(signature, "CTNNB1", "ZNRF3", "APC", "LGR5", "MYC", "TCF7", "TCF7L1", "TCF7L2")
list_signatures <- list(CTNNB1=c("CTNNB1"), ZNRF3=c("ZNRF3"),APC=c("APC"),
                        AXIN2=c("AXIN2"), LGR5=c("LGR5"), SP5=c("SP5"),# other canonical targets of Wnt/B-catenin pathway
                        MYC=c("MYC"),AFF3=c("AFF3"),PDE2A=c("PDE2A"), # other non canonical targets of Wnt/B-catenin pathway
                        LEF1=c("LEF1"), TCF7=c("TCF7"), TCF7L1=c("TCF7L1"),TCF7L2=c("TCF7L2"), # LEF1 and other members of TCF/LEF complex
                        signature=signature, 
                        Nusse_PDE2A=Nusse_PDE2A,
                        KeggWntPathwayGenes=KeggWntPathwayGenes, 
                        pancancer3gene=c("AXIN2","SP5","TNFRSF19") 
) 

# Read inputs: counts (average if multiple probes) and metadata ----------------------------

# Read counts (average if multiple probes) and metadata -----------
Assie <- readRDS(file="./ACC_datasets_formatted/Assie/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")
Assie_metadata <- read.csv("./ACC_datasets_formatted/Assie/Metadata/Metadata_Assie.csv")
Assie_metadata <- Assie_metadata[!(is.na(Assie_metadata$BcatStatus)),]
Assie_metadata$BcatStatus <- factor(Assie_metadata$BcatStatus, levels=c("CTNNB1","ZNRF3","wt"))

Heaton <- readRDS(file="./ACC_datasets_formatted/Heaton/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")
Heaton_metadata <- read.csv("./ACC_datasets_formatted/Heaton/Metadata/Metadata_Heaton.csv")
Heaton_metadata <- Heaton_metadata[Heaton_metadata$Histotype %in% c("ACC"),]
Heaton_metadata <- Heaton_metadata[!(is.na(Heaton_metadata$BetaCateninStaining)),]
Heaton_metadata$BcatStatus <- factor(Heaton_metadata$BetaCateninStaining, levels=c("Nuclear","Membrane"))

tcga <- readRDS(file="./ACC_datasets_formatted/tcga/Log2Pseudocount1_DESeqNormCounts_tcga.RDS")
tcga_Bcat <- readRDS(file="./ACC_datasets_formatted/tcga/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
colnames(tcga_Bcat)[1] <- "TumorID"
tcga_metadata <- read.csv("./ACC_datasets_formatted/tcga/Metadata/Metadata_tcga.csv")
tcga_metadata <- merge(tcga_metadata, tcga_Bcat, all=F, by="TumorID")
tcga_metadata$BcatStatus <- factor(tcga_metadata$Alteration, levels=c("CTNNB1_mut","ZNRF3_mut","wt"))

Lefevre <- readRDS(file="./ACC_datasets_formatted/Lefevre/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")
Lefevre_metadata <- read.csv("./ACC_datasets_formatted/Lefevre/Metadata/Metadata_Lefevre.csv")
Lefevre_metadata <- Lefevre_metadata[Lefevre_metadata$BcatStatus %in% c("expressed","repressed"),] # do not use controls in plots
Lefevre_metadata$BcatStatus <- factor(Lefevre_metadata$BcatStatus, levels=c("expressed","repressed")) 
Lefevre <- Lefevre[,c(1,which(colnames(Lefevre) %in% Lefevre_metadata$Assay.Name) )] # do not use controls in plots

Assie_Bcat <- readRDS(file="./ACC_datasets_formatted/Assie/Metadata/Assie_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
Heaton_Bcat <- readRDS(file="./ACC_datasets_formatted/Heaton/Metadata/Heaton_sampleIDs_BcatNuclearMembrane_onlyACCs.RDS")
tcga_Bcat <- readRDS(file="./ACC_datasets_formatted/TCGA/Metadata/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
Lefevre_Bcat <- Lefevre_metadata[,c('Assay.Name','BcatStatus')]
colnames(Lefevre_Bcat) <- c("Sample","Alteration")

list_datasets <- list(tcga=tcga,Heaton=Heaton,Assie=Assie)

metric <- "sens_constrain"
#metric <- "sum_sens_spec"

# PNG plots with AUC & co values in title ---------------------------------
for(d in 1:length(list_datasets)){
  dataset_name <- names(list_datasets)[d]
  dataset <- list_datasets[[d]]

  print(dataset_name)

  if(dataset_name == "Assie"){
    Bcat_obj <- Assie_Bcat
    alterations <- c("CTNNB1","wt")}
  if(dataset_name == "Heaton"){
    Bcat_obj <- Heaton_Bcat
    alterations <- c("Nuclear","Membrane")}
  if(dataset_name == "tcga"){
    Bcat_obj <- tcga_Bcat
    alterations <- c("CTNNB1_mut","wt")}
  if(dataset_name == "Lefevre"){
    Bcat_obj <- Lefevre_Bcat
    alterations <- c("expressed","repressed")}

  for(s in 1:length(list_signatures)){
    signatureName <- names(list_signatures)[s]
    signature_genes <- list_signatures[[s]]
    print(signatureName)

    ### AUC, sensitivity, specificity values
    if(length(signature_genes)>1){
      signature <- dataset[which(dataset$GeneSymbol %in% signature_genes),-1]
      mean <- apply(as.matrix(signature), 2, mean)
      mean <- as.data.frame(cbind(names(mean),mean))
      colnames(mean) <- c("Sample","Signature")
      input <- merge(mean, Bcat_obj, by="Sample")
      input <- input[input$Alteration %in% alterations,]
      input$Signature <- as.numeric(as.vector(input$Signature))
    }else{
      singleGene <- dataset[which(dataset$GeneSymbol==signature_genes),-1]
      singleGene <- t(rbind(names(singleGene), singleGene))
      colnames(singleGene) <- c("Sample","Signature")
      input <- merge(singleGene, Bcat_obj, by="Sample")
      input <- input[input$Alteration %in% alterations,]
      input$Signature <- as.numeric(as.vector(input$Signature))
    }
    
    if(metric == "sum_sens_spec"){
      cp <- cutpointr(input, Signature, Alteration,
                    method = maximize_metric, metric = sum_sens_spec)
    }
    if(metric == "sens_constrain"){
      cp <- cutpointr(input, Signature, Alteration,
                      method = maximize_metric, metric = sens_constrain)
    }
    #saveRDS(file=paste0(out,"/Cutpointr_", metric, "_", dataset_name, "_", signatureName, ".RDS"), cp)
    
    ### Plot
    plot_out <- plot_Density_byBcatStatus(dataset=dataset, dataset_name=dataset_name, signature=signature_genes, signatureName=signatureName, Bcat=Bcat_obj)
    png(paste0(out,"/GeneSignaturesDensityDistribution_perCTNNB1Status_",dataset_name,"_",signatureName,"_AvgProbes_", metric,".png"))
    print(plot_out +
      ggtitle(paste0(dataset_name, ", AUC=",round(cp$AUC, digits=2), "\nsensitivity=", round(cp$sensitivity,digits=2), ",specificity=", round(cp$specificity,digits=2))))
    dev.off()
  }
}

# HEATMAP AUC ROC values of genes/signatures of interest -------------------------------------------------------
heatmapAUC <- as.data.frame(matrix(ncol=3))
colnames(heatmapAUC) <- c("dataset","signature","AUC")
for(dataset_name in c("Assie","Heaton","tcga")){
  
  for(signatureName in c("CTNNB1","ZNRF3","APC","AXIN2","LGR5","SP5","TCF7","TCF7L1","TCF7L2","KeggWntPathwayGenes", "Nusse_PDE2A", "LEF1", "signature","pancancer3gene")){
    cp <- readRDS(file=paste0(out,"/Cutpointr_", metric, "_", dataset_name, "_", signatureName, ".RDS"))
    vec <- c(dataset_name, signatureName, round(cp$AUC, digits=3)) 
    names(vec) <- c("dataset","signature","AUC")
    heatmapAUC <- rbind(heatmapAUC, vec)
  }
  
}
heatmapAUC <- heatmapAUC[-1,]
heatmapAUC$AUC <- as.numeric(as.vector(heatmapAUC$AUC))
saveRDS(file=paste0(out, "HeatmapAUCs_dataset_signature_",metric,".RDS"), heatmapAUC)

heatmapAUC <- readRDS(file=paste0(out, "HeatmapAUCs_dataset_signature_",metric,".RDS"))
heatmapAUC <- heatmapAUC[heatmapAUC$signature %in% c("CTNNB1","ZNRF3","APC","AXIN2","LGR5","SP5","TCF7","TCF7L1","TCF7L2","KeggWntPathwayGenes", "Nusse_PDE2A", "LEF1", "signature"),]
heatmapAUC$signature <- factor(heatmapAUC$signature, levels=c("CTNNB1","ZNRF3","APC","AXIN2","LGR5","SP5","TCF7","TCF7L1","TCF7L2","LEF1","KeggWntPathwayGenes", "Nusse_PDE2A", "signature"))
png(paste0(out, "HeatmapAUCs_dataset_GenesSignature_ForSuppFig_",metric, ".png"))
ggplot(data=heatmapAUC, aes(x=dataset,y=signature,fill=AUC)) + geom_tile() +
  geom_text(aes(label=AUC)) + theme_bw() +
  scale_fill_gradient2(low="navy", mid="white", high="red", midpoint = 0.6)
dev.off()
