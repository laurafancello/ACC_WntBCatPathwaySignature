library(GEOquery)
library(ggplot2)
library(ggrepel)
source("./Scripts/Miscellaneous_Functions.R")
outPath <- "./DEanalyses_ACC/CTNNB1mutVsWt/"

#  IMPORT DATASET DIRECTLY FROM GEO ---------------------------------------
# load series and platform data from GEO
gset <- getGEO("GSE10927", GSEMatrix =TRUE, AnnotGPL=TRUE)
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
          exprs(gset) <- log2(ex) 
          ex <- log2(ex)}

# Read microarray annotation (probeID - geneID) ---------------------------
microarray_annotation <- read.table(file="./ACC_datasets_formatted/MicroarrayAnnotations/GPL570.an.txt", sep="\t", header=T, fill=T, quote="") # probeID to geneSymbol file from https://gemma.msl.ubc.ca/arrays/showArrayDesign.html?id=4
colnames(microarray_annotation)[c(1,2)] <- c("FROM","TO")

ex_out <- as.data.frame(cbind(rownames(ex), ex))
colnames(ex_out)[1] <- 'FROM'
ex_out <- merge(ex_out, microarray_annotation[,c('FROM','TO')], by='FROM', all=F)
ex_out <- ex_out[,c(1,ncol(ex_out),2:(ncol(ex_out)-1))]
colnames(ex_out)[1:2] <- c('Probe','Gene')

#  Read metadata and classify samples -------------------------------------
metadata <- read.csv(file="./ACC_datasets_formatted/Heaton/Metadata/Metadata_Heaton.csv")
metadata <- metadata[metadata$Histotype %in% c("ACC"),] # Keep only ACCs in metadata object

# Add a column specifying the condition of the sample, by which to perform DE analysis
metadata <- metadata[!(is.na(metadata$BetaCateninStaining)),]
metadata$condition <- NA
metadata[metadata$BetaCateninStaining %in% c("Nuclear"),]$condition <- "altered"
metadata[metadata$BetaCateninStaining %in% c("Membrane"),]$condition <- "wt"
# Factorize the newly added column
metadata$condition <- factor(metadata$condition, levels=c("wt","altered"))

# Keep only ACCs in counts object
ex_filt <- ex[,which(colnames(ex) %in% metadata$GSM_ID)]

# Check that counts and metadata samples order is the same
ex_filt_ok <- ex_filt[,match(metadata$GSM_ID,colnames(ex_filt))]

if(identical(colnames(ex_filt_ok), metadata$GSM_ID )){

  # Perform DE analysis -----------------------------------------------------
  DE_result <- DE_microarray(expr=ex_filt_ok, metadata=metadata, condition=condition, microarray_annotation=microarray_annotation)
  res <- as.data.frame(DE_result$affyDEResults)
  write.csv(file=paste0(outPath, "DEanalysis_CTNNB1onlyMut_vs_Wt_allGenes_Heaton_ACCs.csv"), res)
}  
