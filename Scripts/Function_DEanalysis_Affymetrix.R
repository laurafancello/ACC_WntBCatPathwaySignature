################################################################################
# Define function for DE analysis with microarray data ------------------------------
################################################################################
DE_microarray <- function(expr, metadata, condition, microarray_annotation){
  levelsCondition <- c("wt","altered")
  if(!(identical(levels(metadata$condition), levelsCondition))){
    print("condition must be a factor of metadata object with only two levels: wt and altered")
  }
  # Create the SummarizedExperiment object
  affyDataset <- SummarizedExperiment::SummarizedExperiment(assays = expr, colData = metadata)
  # Create a design matrix
  affyDesign <- model.matrix(~0 + condition, data = SummarizedExperiment::colData(affyDataset))
  # Avoid special characters in column names
  colnames(affyDesign) <- make.names(colnames(affyDesign))
  # Create a constrast matrix
  affyContrast <- limma::makeContrasts(conditionaltered-conditionwt, levels=affyDesign)
  # Run differential expression analysis
  affyDEExperiment <- RCPA::runDEAnalysis(affyDataset, method = "limma", design = affyDesign, contrast = affyContrast, annotation = microarray_annotation)
  # Extract the differential analysis result
  affyDEResults <- SummarizedExperiment::rowData(affyDEExperiment)
  return(list(affyDataset=affyDataset,affyDesign=affyDesign,affyContrast=affyContrast,affyDEExperiment=affyDEExperiment,affyDEResults=affyDEResults))
}