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


# SURVIVAL ANALYSIS TWO GROUPS -------------------------------------------------------
print("------ Start survival analyses -----------------------------------------")
for(d in 1:length(list_datasets)){
  dataset_name <- names(list_datasets)[d]
  dataset <- list_datasets[[d]]

  print(dataset_name)

  clin <- read.csv(file=paste0("C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/", dataset_name, "/Metadata/Metadata_", dataset_name, ".csv"))
  clin <- clin[clin$Histotype == "ACC",] # select only ACC samples (required only for Heaton)

  for(s in 1:length(list_signatures)){
    signature_name <- names(list_signatures)[s]
    signature_genes <- list_signatures[[s]]
    print(signature_name)

      for(threshold in c("cutpointr","median")){

        print(threshold)

        if(dataset_name == "tcga"){
          log2counts <- as.data.frame(readRDS(file=paste0("C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/", dataset_name, "/Log2Pseudocount1_DESeqNormCounts_tcga.RDS")))
        }else{
          log2counts <- as.data.frame(readRDS(file=paste0("C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/", dataset_name, "/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")))
          GeneSymbol <- log2counts$GeneSymbol
          log2counts <- log2counts[,which(colnames(log2counts) %in% clin$TumorID)]
          log2counts <- as.data.frame(cbind(GeneSymbol, log2counts))
        }

        signature <- log2counts[log2counts$GeneSymbol %in% signature_genes,-1]
        if(length(signature_genes)>1){
          mean <- apply(as.matrix(signature), 2, mean)
          mean <- as.data.frame(cbind(names(mean),mean))
        }else{
          mean <- as.data.frame(t(rbind(colnames(signature),signature)))
        }
        colnames(mean) <- c("patientId","Signature")
        mean$Signature <- as.numeric(as.vector(mean$Signature))

        if(threshold=="cutpointr"){
          cp <- readRDS(file=paste0(out, "Cutpointr_sens_constrain_bootstrap_", n_boot, "_", dataset_name, "_", signature_name,".RDS"))
          cutoff <- round(cp$optimal_cutpoint, digits=2)
        }
        if(threshold=="median"){
          cutoff <- round(median(mean$Signature), digits=2)
        }

        low <- mean[mean$Signature<cutoff,]
        high <- mean[mean$Signature>=cutoff,]

        low <- as.data.frame(cbind(low, rep("low", length(low$patientId))))
        colnames(low)[ncol(low)] <- "SignatureLevel"
        high <- as.data.frame(cbind(high, rep("high", length(high$patientId))))
        colnames(high)[ncol(high)] <- "SignatureLevel"
        SignatureLevels <- rbind(low, high)

        SignatureLevels$TumorID <- SignatureLevels$patientId

        input_surv <- merge(SignatureLevels, clin, by="TumorID", all=F)
        if(length(which(input_surv$OS_Status %in% c("alive", "dead")))>1){
          input_surv[input_surv$OS_Status %in% c("alive"),]$OS_Status <- 0
          input_surv[input_surv$OS_Status %in% c("dead"),]$OS_Status <- 1
        }
        input_surv$OS_Status <- as.numeric(as.vector(input_surv$OS_Status))

        input_surv$SignatureLevel <- factor(input_surv$SignatureLevel, levels=c("low","high"))

        fit <- survfit(Surv(input_surv$OS_Months, input_surv$OS_Status) ~
                         input_surv$SignatureLevel, data = input_surv)
        values <- surv_pvalue(fit)
        p <- round(values$pval, digits=3)

        # PNG
        png(file = paste0(out, "Survival/", dataset_name, "_survival_curve_", signature_name, "_",threshold,"_publi.png"), width=600, height = 600)
        print(ggsurvplot(fit=fit,
                         linetype = c(1,3), censor.shape = 124, censor.size=2,
                         legend.title = "", font.x = 22, font.y = 22, title = paste0(dataset_name, ", ", signature_name, ", ", threshold, "\n p=", p),
                         font.tickslab = 18, font.legend = 18, pval.size = 7,
                         legend.labs = c("high","low"),
                         palette = c("red", "blue"), surv.scale = "percent",
                         ylab = "% of survival", xlab = "Time (in months)",
                         ggtheme = theme_juju(), risk.table = T, pval = T,
                         tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                         risk.table.y.text.col = F))
        dev.off()
        # PDF
        pdf(file = paste0(out, "Survival/", dataset_name, "_survival_curve_", signature_name, "_",threshold,"_publi.pdf"), useDingbats = F)
        print(ggsurvplot(fit=fit,
                         linetype = c(1,3), censor.shape = 124, censor.size=2,
                         legend.title = "", font.x = 20, font.y = 20, title =  paste0(dataset_name, ", ", signature_name, ", ", threshold, "\n p=", p),
                         font.tickslab = 18, font.legend = 18, pval.size = 7,
                         legend.labs = c("high","low"),
                         palette = c("red", "blue"), surv.scale = "percent",
                         ylab = "% of survival", xlab = "Time (in months)",
                         ggtheme = theme_juju(), risk.table = T, pval = T,
                         tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                         risk.table.y.text.col = F))
        dev.off()
      }
    }
}


# SURVIVAL ANALYSIS THREE GROUPS -------------------------------------------------------
print("------ Start survival analyses -----------------------------------------")
for(d in 1:length(list_datasets)){
  dataset_name <- names(list_datasets)[d]
  dataset <- list_datasets[[d]]
  
  print(dataset_name)
  
  clin <- read.csv(file=paste0("C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/", dataset_name, "/Metadata/Metadata_", dataset_name, ".csv"))
  clin <- clin[clin$Histotype == "ACC",] # select only ACC samples (required only for Heaton)
  
  for(s in 1:length(list_signatures)){
    signature_name <- names(list_signatures)[s]
    signature_genes <- list_signatures[[s]]
    print(signature_name)
      
    if(dataset_name == "tcga"){
      log2counts <- as.data.frame(readRDS(file=paste0("C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/", dataset_name, "/Log2Pseudocount1_DESeqNormCounts_tcga.RDS")))
    }else{
      log2counts <- as.data.frame(readRDS(file=paste0("C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/", dataset_name, "/Log2NormData_all_SampleID_ok_AvgMultipleProbesAnyGene.RDS")))
      GeneSymbol <- log2counts$GeneSymbol
      log2counts <- log2counts[,which(colnames(log2counts) %in% clin$TumorID)]
      log2counts <- as.data.frame(cbind(GeneSymbol, log2counts))
    }
    
    signature <- log2counts[log2counts$GeneSymbol %in% signature_genes,-1]
    if(length(signature_genes)>1){
      mean <- apply(as.matrix(signature), 2, mean)
      mean <- as.data.frame(cbind(names(mean),mean))
    }else{
      mean <- as.data.frame(t(rbind(colnames(signature),signature)))
    }
    colnames(mean) <- c("patientId","Signature")
    mean$Signature <- as.numeric(as.vector(mean$Signature))
    
    quantiles <- quantile(mean$Signature, probs = c(0,0.33,0.66,1))
    
    low <- mean[mean$Signature<=quantiles[2],]
    mid <- mean[((mean$Signature>quantiles[2])&(mean$Signature<quantiles[3])),]
    high <- mean[mean$Signature>=quantiles[3],]
    
    low <- as.data.frame(cbind(low, rep("low", length(low$patientId))))
    colnames(low)[ncol(low)] <- "SignatureLevel"
    mid <- as.data.frame(cbind(mid, rep("mid", length(mid$patientId))))
    colnames(mid)[ncol(mid)] <- "SignatureLevel"
    high <- as.data.frame(cbind(high, rep("high", length(high$patientId))))
    colnames(high)[ncol(high)] <- "SignatureLevel"
    SignatureLevels <- rbind(low,mid,high)
    
    SignatureLevels$TumorID <- SignatureLevels$patientId
    
    input_surv <- merge(SignatureLevels, clin, by="TumorID", all=F)
    if(length(which(input_surv$OS_Status %in% c("alive", "dead")))>1){
      input_surv[input_surv$OS_Status %in% c("alive"),]$OS_Status <- 0
      input_surv[input_surv$OS_Status %in% c("dead"),]$OS_Status <- 1
    }
    input_surv$OS_Status <- as.numeric(as.vector(input_surv$OS_Status))
    
    input_surv$SignatureLevel <- factor(input_surv$SignatureLevel, levels=c("low","mid","high"))
    
    fit <- survfit(Surv(input_surv$OS_Months, input_surv$OS_Status) ~
                     input_surv$SignatureLevel, data = input_surv)
    values <- surv_pvalue(fit)
    p <- round(values$pval, digits=3)
    
    # PNG
    png(file = paste0(out, "Survival/", dataset_name, "_survival_curve_", signature_name, "_3groups_33_66Quantiles_publi.png"), width=600, height = 600)
    print(ggsurvplot(fit=fit,
                     #linetype = c(1,3,1), 
                     censor.shape = 124, censor.size=2,
                     legend.title = "", font.x = 22, font.y = 22, title = paste0(dataset_name, ", ", signature_name, ", 3groups_33_66Quantiles \n p=", p),
                     font.tickslab = 18, font.legend = 18, pval.size = 7,
                     #legend.labs = c("high","mid","low"),
                     #palette = c("red", "gray30","blue"), 
                     surv.scale = "percent",
                     ylab = "% of survival", xlab = "Time (in months)",
                     ggtheme = theme_juju(), risk.table = T, pval = T,
                     tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                     risk.table.y.text.col = F))
    dev.off()
    # PDF
    pdf(file = paste0(out, "Survival/", dataset_name, "_survival_curve_", signature_name, "_3groups_33_66Quantiles_publi.pdf"), useDingbats = F)
    print(ggsurvplot(fit=fit,
                     #linetype = c(1,3,1), 
                     censor.shape = 124, censor.size=2,
                     legend.title = "", font.x = 22, font.y = 22, title = paste0(dataset_name, ", ", signature_name, ", 3groups_33_66Quantiles \n p=", p),
                     font.tickslab = 18, font.legend = 18, pval.size = 7,
                     #legend.labs = c("high","mid","low"),
                     #palette = c("red", "gray30","blue"), 
                     surv.scale = "percent",
                     ylab = "% of survival", xlab = "Time (in months)",
                     ggtheme = theme_juju(), risk.table = T, pval = T,
                     tables.theme = theme_cleantable_juju(), fontsize = 5.5,
                     risk.table.y.text.col = F))
    dev.off()
  }
}


