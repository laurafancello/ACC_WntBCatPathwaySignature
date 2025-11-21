library(dplyr)
library(stringr)

# Counts ------------------------------------------------------------------
norm <- readRDS("C:/Users/LF260934/Documents/WORKING_FOLDER/Results/DEanalyses_ACC_LEF1highLow_multipleDatasets/tcga/TCGA_DESeqNormalizedCounts.RDS")
norm_log2 <- as.data.frame(apply(norm, 2, function(x) log2(x+1)))
norm_log2 <- as.data.frame(cbind(rownames(norm_log2), norm_log2))
colnames(norm_log2)[1] <- "GeneSymbol"
saveRDS(file="C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/tcga/Log2Pseudocount1_DESeqNormCounts_tcga.RDS", norm_log2)
 


#  Metadata ---------------------------------------------------------------
BcatStatus <- readRDS(file="C:/Users/LF260934/Documents/DATASETS/TCGA/tcga_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")
colnames(BcatStatus) <- c("TumorID","BcatStatus")
#TP53Status <- readRDS(file="C:/Users/LF260934/Documents/DATASETS/TCGA/tcga_sampleIDs_TP53altered.RDS")
TP53Status <- read.table(file="C:/Users/LF260934/Documents/WORKING_FOLDER/Results/TCGA_ACC/acc_tcga_study_17062024/ACC_TCGA_mut_CN_TP53.txt", header=T, sep="\t")
colnames(TP53Status)[c(1,ncol(TP53Status))] <- c("TumorID","TP53Status")
clin <- read.csv(file="C:/Users/LF260934/Documents/WORKING_FOLDER/Results/TCGA_ACC/acc_tcga_study_17062024/TCGA_acc_cli.csv")
clin <- clin[,-1]
colnames(clin)[2] <- "TumorID"
all <- merge(BcatStatus, clin, by="TumorID", all=T)
all  <- merge(all, TP53Status, by="TumorID", all=T)
all <- all %>% dplyr::mutate(OS5years = ifelse(OS_MONTHS <= 60 & OS_STATUS == "1:DECEASED", 1, 0))
#all[is.na(all$OS_STATUS),]$OS5years <- NA
all$Histotype <- "ACC"
index_first <- which(colnames(all) %in% c("TumorID","Histotype","BcatStatus","patientId","AGE","DFS_MONTHS", "DFS_STATUS", "OS_MONTHS", "OS_STATUS", 
                                          "HISTORY_ADRENAL_HORMONE_EXCESS",
                                          "HISTOLOGICAL_DIAGNOSIS","SEX","WEISS_SCORE_OVERALL", "TP53Status", "OS5years", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                                          "SAMPLE_INITIAL_WEIGHT"))                               
index_last <- which(!(colnames(all) %in% c("TumorID","Histotype","BcatStatus","patientId","AGE","DFS_MONTHS", "DFS_STATUS", "OS_MONTHS", "OS_STATUS", 
                                           "HISTORY_ADRENAL_HORMONE_EXCESS",
                                          "HISTOLOGICAL_DIAGNOSIS","SEX","WEISS_SCORE_OVERALL", "TP53Status", "OS5years", "AJCC_PATHOLOGIC_TUMOR_STAGE",
                                          "SAMPLE_INITIAL_WEIGHT")))                           
all_ok <- all[,c(index_first, index_last)]
all_ok$HormoneExpression <- all$HISTORY_ADRENAL_HORMONE_EXCESS
all_ok$HormoneExpression <- stringr::str_replace_all(all_ok$HormoneExpression, "\\|", "+")
all_ok$HormoneExpression <- stringr::str_replace_all(all_ok$HormoneExpression, "Androgen+Cortisol", "Cortisol+Androgen")

### Format like other datasets overall survival data
all_ok$OS_Status <- stringr::str_replace_all(all_ok$OS_STATUS,"1:DECEASED","1")
all_ok$OS_Status <- stringr::str_replace_all(all_ok$OS_Status,"0:LIVING","0")
all_ok$OS_Status <- as.numeric(as.vector(all_ok$OS_Status))
all_ok$OS_Months <- as.numeric(as.vector(all_ok$OS_MONTHS))

write.csv(file="C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/tcga/Metadata/Metadata_tcga.csv", all_ok, row.names = F)

