inPath <- "C:/Users/LF260934/Documents/DATASETS/ACC_AssieNatGenet2014_microarray/"

Supp2A <- read.table(file=paste0(inPath,"/Metadata/Supp2A.txt"), sep="\t", header=T)
Supp3 <- read.table(file=paste0(inPath,"/Metadata/Supp3.txt"), sep="\t", header=T)
Supp5A <- read.table(file=paste0(inPath,"/Metadata/Supp5A.txt"), sep="\t", header=T)
Supp5B <- read.table(file=paste0(inPath,"/Metadata/Supp5B.txt"), sep="\t", header=T)
Supp1B <- read.csv(file=paste0(inPath,"/Metadata/SuppTable1b.csv"))
colnames(Supp1B)[2] <- "TumorID"

metadata <- merge(Supp2A, Supp3, by="TumorID", all=T)
metadata <- merge(metadata, Supp5A[,c("TumorID","mRNAsubgroup","DNAmethylationSubgroup","miRNAcluster","mir.506.514","DLK1.MEG3")], by="TumorID", all=T)
metadata <- merge(metadata, Supp5B, by="TumorID", all=T)

metadata <- merge(metadata, Supp1B[,c("TumorID","Sex","Age.at.diagnosis","Tumor.size..mm.","Hormonal.Secretion")], by="TumorID", all=T)
fileTumor <- read.table(file=paste0(inPath,"/Metadata/FileID_TumorID.txt"), sep="\t", header=F)
colnames(fileTumor) <- c("GSM_ID","TumorID")
metadata_microarray <- merge(metadata, fileTumor, by="TumorID", all.x=F, all.y=T)

sampleId_classification <- readRDS(file="C:/Users/LF260934/Documents/DATASETS/ACC_AssieNatGenet2014_microarray/Metadata/Assie_sampleIDs_CTNNB1mut_ZNRF3mut_wt.RDS")


# Format metadata with field names according to my common nomenclature
metadata_formatOK <- metadata_microarray
colnames(metadata_formatOK) <- c("Sample","NbMutations","OS5years","WeissScore","ENSATstage","ZNRF3","CTNNB1","TP53","CDKN2A","RB1","MEN1","DAXX","MED12",
                                 "TERT","Wnt.pathway","p53.pathway","Chromatin.pathway","mRNAsubgroup","DNAmethylationSubgroup","miRNAcluster","mir.506.514","DLK1.MEG3",
                                 "OS_Status","OS_Months","MolecularGroup", "Sex","Age","Tumor.size..mm.","HormoneExpression","GSM_ID")
metadata_formatOK <- merge(metadata_formatOK, sampleId_classification, by="Sample", all=T)
colnames(metadata_formatOK)[ncol(metadata_formatOK)] <- "CTNNB1_ZNRF3_Alteration" # change name of the last column, which was introduce by merge at previous line

# Add microarray platform (GPLXXX) 
metadata_formatOK$MicroarrayPlatform <- "GPL16686"

# Add histotype (even if ACC for all)
metadata_formatOK$Histotype <- "ACC"
metadata_formatOK$TumorID <- metadata_formatOK$Sample

metadata_formatOK$HormoneExpression <- stringr::str_replace(metadata_formatOK$HormoneExpression, "no", "None")

# Save
saveRDS(file="C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/Assie/Metadata/Metadata_Assie.RDS", metadata_formatOK)
write.csv(file="C:/Users/LF260934/Documents/WORKING_FOLDER/ACC_datasets_formatted/Assie/Metadata/Metadata_Assie.csv", metadata_formatOK, row.names = F)
