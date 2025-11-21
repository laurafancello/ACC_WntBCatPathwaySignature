jouinot=readRDS(file="C:/Users/LF260934/Documents/DATASETS/ACC_Cochin_bulkFFPE/RNAseq3p_bulk_FFPE_EJE2022_log2DESeqPseudocount1.RDS")
jouinot <- as.data.frame(jouinot)
jouinot <- as.data.frame(cbind(rownames(jouinot), jouinot))
colnames(jouinot)[1] <- "GeneSymbol"
saveRDS(file="C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/JouinotFFPE/Log2NormData_all_SampleID_ok.RDS",jouinot)
