library(GEOquery)
library(ggplot2)
library(ggrepel)
source("./Scripts/Function_DEanalysis_Affymetrix.R")
source("./Scripts/DEanalysisAndCo_Functions.R")
inPath <- "./ACC_datasets_formatted/"
outPath <- "./DEanalyses_ACC/CTNNB1mutVsWt/"

#  IMPORT DATASET DIRECTLY FROM GEO ---------------------------------------
# load series and platform data from GEO
gset <- getGEO("GSE49278", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))


# LOG2 TRANSFORMATION -----------------------------------------------------
# log2 transformation, if necessary 
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN # Assie data were already log2 transformed, so no need to
  ex <- log2(ex)
}

# Read the metadata files (SOURCE DATA FROM FIGURES IN ASSIE et al) --------
metadata <- read.csv(file=paste0(inPath,"/Assie/Metadata/Metadata_Assie.csv"))

# PREPARE AND PERFORM DE ANALYSIS CTNNB1mut VS WT -------------------------
norm_data_ok <- ex
# Format to be compatible with DE-microarray function
metadata <- metadata[metadata$CTNNB1_ZNRF3_Alteration %in% c("CTNNB1","wt"),] # keep only CTNNB1 and wt in metadata
norm_data_ok <- norm_data_ok[,which(colnames(norm_data_ok) %in% metadata$GSM_ID)] # keep only CTNNB1 and wt in counts
metadata$condition <- ifelse(metadata$CTNNB1_ZNRF3_Alteration == "CTNNB1", "altered", "wt")
metadata$condition <- factor(metadata$condition, levels=c("wt","altered"))

# Check that counts and sample_info objects have same sample order
norm_data_ok <- norm_data_ok[,match(metadata$GSM_ID, colnames(norm_data_ok))] 
identical(colnames(norm_data_ok), metadata$GSM_ID)

# Read microarray annotation (probeID - geneID)
microarray_annotation <- read.table(file=paste0(inPath,"/MicroarrayAnnotations/GPL16686.an.txt"), sep="\t", header=T, fill=T, quote="") # probeID to geneSymbol file from https://gemma.msl.ubc.ca/arrays/showArrayDesign.html?id=787
colnames(microarray_annotation)[c(1,2)] <- c("FROM","TO")

DE_result <- DE_microarray(expr=norm_data_ok, metadata=metadata, condition=condition, microarray_annotation=microarray_annotation)
res <- DE_result$affyDEResults
write.csv(file=paste0(outPath, "DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_Assie.csv"), res) # Save results

# DE analysis without ACC12,ACC13,ACC18,ACC36,ACC38 (CTNNB1wt but LEF1 and signature high) --------
#Redo differential expression analysis CTNNB1 mut vs wt after removal of  ACC12,ACC13,ACC18,ACC36,ACC38 which are CTNNB1_wt but
# LEF1/B-catenin signature high to see if I am able to find significant DE B-catenin targets
remove <- metadata[metadata$TumorID %in% c("ACC12","ACC13","ACC18","ACC36","ACC38"),]$GSM_ID
metadata_rm <- metadata[-(which(metadata$GSM_ID %in% remove)),]
norm_data_ok_rm <- norm_data_ok[-which(colnames(norm_data_ok) %in% remove),]
norm_data_ok_DEF <- norm_data_ok_rm[,match(metadata_rm$GSM_ID, colnames(norm_data_ok_rm))] 
identical(colnames(norm_data_ok_DEF), metadata_rm$GSM_ID)

DE_result <- DE_microarray(expr=norm_data_ok_DEF, metadata=metadata_rm, condition=condition, microarray_annotation=microarray_annotation)
res <- as.data.frame(DE_result$affyDEResults)
write.csv(file=paste0(outPath, "DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_Assie_withoutACCmisclassified.csv"), res)



# COMPARE FDR/FC OF B-CATENIN TARGETS IN THE 2 DE ANALYSIS ----------------
Nusse_PDE2A <- c("ABCB1", "AFF3", "AXIN2", "BCL2L2", "BIRC5", "CCND1", "CDC25A", "CDKN2A", "CDX1", "CLDN1", "CTLA4", "DKK1", "EDN1", "ENAH", "ENC1", "FGF18", 
                 "FGF4", "FGFBP1", "FOSL1", "FSCN1", "FST", "FZD7", "GBX2", "HES1", "HNF1A", "ID2", "JAG1", "JUN", "KRT5", "L1CAM","LAMC2", "LEF1", "LGR5", "MMP14", 
                 "MMP7", "MYC", "MYCBP", "NEDD9", "NEUROD1", "NEUROG1", "NOS2", "NOTCH2", "NRCAM", "PDE2A", "PLAU", 
                 "PPARD", "PTGS2", "S100A4", "SGK1", "SMC3", "SP5", "SUZ12", "TBX3", "TBXT", "TCF4", "TERT", "TNC", "TNFRSF19", "VCAN", "VEGFA", 
                 "YY1AP1") # known B-catenin targets (those reported on Nusse lab website + Herbst et al 2014 paper + PDE2A + AFF3)
ACC_BcatSignature <- c("AFF3","AXIN2","FSCN1","LEF1","MYC","PDE2A","SP5","TBX3","TNFRSF19")

res <- read.csv(file=paste0(outPath, "DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_Assie.csv"))
res_noMiscalssified <- read.csv(file=paste0(outPath, "DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_Assie_withoutACCmisclassified.csv"))

sub_res <- res[res$ID %in% Nusse_PDE2A, c("ID","logFC","pFDR")]
colnames(sub_res) <- c("Gene","logFC_allTumors","pFDR_allTumors")
sub_res_noMisclassified <- res_noMiscalssified[res_noMiscalssified$ID %in% Nusse_PDE2A, c("ID","logFC","pFDR")]
colnames(sub_res_noMisclassified) <- c("Gene","logFC_noMisclassified","pFDR_noMisclassified")

merged <- merge(sub_res, sub_res_noMisclassified, by="Gene")
merged$In_signature <- ifelse(merged$Gene %in% ACC_BcatSignature, "TRUE", "FALSE")

merged <- merged[order(merged$pFDR_noMisclassified),]
write.csv(file=paste0(outPath, "Compare_FC_FDR_BCateninTargets_WithOrWithoutACCmisclassified.csv"), merged, row.names = F)

