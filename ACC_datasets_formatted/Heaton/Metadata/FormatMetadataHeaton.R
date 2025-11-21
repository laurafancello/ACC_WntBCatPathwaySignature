library(stringr)
metadata <- read.csv(file="C:/Users/LF260934/Documents/DATASETS/ACC_HeatonAmJPathol2012_microarray/Metadata/Metadata.csv")
metadata <- t(metadata)
metadata <- as.data.frame(cbind(rownames(metadata), metadata))
colnames(metadata) <- metadata[1,]
metadata <- metadata[-1,] # remove row which became header
metadata <- metadata[,c(-4,-7)] # rtemove columns with not useful/redundant info
colnames(metadata) <- c("GSM_ID","TumorID","Histotype","CELfile","SourceName","Age","Sex","ClinicalCharacteristics",
                        "TumorDiameter_cm","TumorWeight_g","WeissScore","MitoticRate","TumorStage","OS_Years","OS_Status",
                        "BetaCateninStaining","CTNNB1mut","HormoneExpression")
metadata <- metadata[-nrow(metadata),] # remove last row which is empty
metadata <- metadata[,c(2:ncol(metadata),1)] # put GSM_ID in last position to respect column order expected by previously written code
metadata$Histotype <- stringr::str_replace_all(metadata$Histotype, "normal", "NORMAL") # replace normal by NORMAL to respect values expcted by previously written code
metadata$CTNNB1mut <- stringr::str_replace_all(metadata$CTNNB1mut, "unknown", "NA") # replace unknown by NA to respect values expcted by previously written code
metadata$BetaCateninStaining <- stringr::str_replace_all(metadata$BetaCateninStaining, "unknown", "NA") # replace unknown by NA to respect values expcted by previously written code
metadata$CTNNB1mut <- stringr::str_replace(metadata$CTNNB1mut, "wild type", "wt") # replace to respect values expcted by previously written code
metadata$CTNNB1mut <- stringr::str_replace(metadata$CTNNB1mut, "mutant", "mut") # replace to respect values expcted by previously written code

metadata$OS_Months <- as.numeric(as.vector(metadata$OS_Years))*12
metadata$MicroarrayPlatform <- "GPL570"
metadata$HormoneExpression <- stringr::str_replace_all(metadata$HormoneExpression, "Aldosterone", "Mineralocorticoids")
metadata$HormoneExpression <- stringr::str_replace_all(metadata$HormoneExpression, "Testosterone", "Androgen")
metadata$HormoneExpression <- stringr::str_replace_all(metadata$HormoneExpression, "No_information", "NA")
metadata$HormoneExpression <- stringr::str_replace_all(metadata$HormoneExpression, "unknown", "NA")
metadata$HormoneExpression <- stringr::str_replace_all(metadata$HormoneExpression, "Non-functional", "None")
metadata$HormoneExpression <- stringr::str_replace_all(metadata$HormoneExpression, "")
write.csv(file="C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/Heaton/Metadata/Metadata_Heaton.csv", metadata, row.names=F)
