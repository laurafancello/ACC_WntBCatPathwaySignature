# Function heatmaps significant GSEA enrichments -----------------------------------
# To compare significant GSEA enrichments between different datasets/patients/cell types etc
# Arguments:
# geneSetsSuffix = suffix of GSEA output file indicating gene set DB used
# outpath = path to output
# grouping1 = first criterium to separate data (heatmap x axis)
# grouping2 = first criterium to separate data (heatmap facet in facet_grid/facet_wrap)

heatmapCommonGSEAEnrichments <- function(geneSetsSuffix, grouping1, grouping2, Correlation_type){
  ### Plot heatmaps
  index1 <- grep(grouping1, colnames(sub_pos))
  index2 <- grep(grouping2, colnames(sub_pos))
  sub_pos <- all_pos[all_pos$GeneSet == geneSetsSuffix,]
  pos <- ggplot(dat=sub_pos, aes(y=pathway, x=sub_pos[,index1])) + geom_tile() + facet_wrap(~sub_pos[,index2]) +
    theme_classic() +theme(axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("dataset") +
    ggtitle(paste0("GSEA ",geneSet," significant (padj<0.1) positive enrichments\non ",Correlation_type,"correlation B-catenin signature- gene expression"))
  sub_neg <- all_neg[all_neg$GeneSet == geneSet,]
  
  neg <- ggplot(dat=sub_neg, aes(y=pathway, x=sub_neg[,index1])) + geom_tile() + facet_wrap(~sub_neg[,index2]) +
    theme_classic() +theme(axis.text.x = element_text(angle=45, hjust=1)) +
    xlab("dataset") +
    ggtitle(paste0("GSEA ",geneSet," significant (padj<0.1) negative enrichments\non ",Correlation_type,"correlation B-catenin signature- gene expression"))
  
  return(list(positiveEnrichmentsHeatmap=pos, negativeEnrichmentsHeatmap=neg))
}




# Density plot mut vs wt by gene/signature expression ---------------------
# Function to plot density of mean signature expression (or gene expression) by Bcat status (CTNNB1mut vs wt or nuclear vs membrane) --------------------
plot_Density_byBcatStatus <- function(dataset, dataset_name, signature, signatureName, Bcat){
  colnames(Bcat) <- c("TumorID","BcatStatus")
  
  select <- dataset[dataset$GeneSymbol %in% signature,]
  select <- as.matrix(select[,-1])
  select <- apply(select, 2, as.numeric)
  if(length(signature)>1){
    mean <- apply(select, 2, mean)
    select <- as.data.frame(cbind(names(mean), mean))
  }else{select <- as.data.frame(cbind(names(select), select))}
  colnames(select) <- c("TumorID","BcatSignature")
  
  all <- merge(select, Bcat, by="TumorID")
  all$BcatSignature <- as.numeric(as.vector(all$BcatSignature))
  all$BcatStatus <- as.factor(all$BcatStatus)
  
  if(length(intersect(unique(all$BcatStatus), c("CTNNB1","ZNRF3","wt")))==3){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1","ZNRF3","wt"))
    colors_manual <- c('red','orange', 'green4')
  }
  if((length(unique(all$BcatStatus))==3)&(length(intersect(unique(all$BcatStatus), c("CTNNB1_mut","ZNRF3_mut","wt")))==3)){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1_mut","ZNRF3_mut","wt"))
    colors_manual <- c('red','orange', 'green4')
  }
  if((length(unique(all$BcatStatus))==4)&(length(intersect(unique(all$BcatStatus), c("CTNNB1_mut","ZNRF3_mut","ZNRF3del","wt")))==4)){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1_mut","ZNRF3_mut","ZNRF3del","wt"))
    colors_manual <- c('red','orange', 'gold1', 'green4')
  }
  if(length(intersect(unique(all$BcatStatus), c("Nuclear","Membrane")))==2){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("Nuclear","Membrane"))
    colors_manual <- c('red', 'green4')
  }
  if(length(intersect(unique(all$BcatStatus), c("CTNNB1_mut","ZNRF3_mut","wt")))==2){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("CTNNB1_mut","wt"))
    colors_manual <- c('red', 'green4')
  }
  if(length(intersect(unique(all$BcatStatus), c("mut","wt")))==2){
    all$BcatStatus <- factor(all$BcatStatus, levels=c("mut","wt"))
    colors_manual <- c('red', 'green4')
  }
  #Plot density distribution B-catenin activation signature separating CTNNB1 mutated vs wt samples
  ggplot(data=all, aes(x=BcatSignature, group=BcatStatus, fill=BcatStatus, color=BcatStatus)) +
    geom_density(alpha=0.5) +
    theme_bw() +
    xlab(signatureName) +
    ggtitle(dataset_name)
}



# Functions survival plots as in Justine paper ----------------------------
theme_juju <- function(){
  theme_classic()+
    theme(plot.title = element_text(hjust=0.5, face = "plain", size = 24),
          axis.text = element_text(color="black"),
          legend.title = element_text(hjust=0.5))
}

theme_cleantable_juju <- function(){
  theme_cleantable()+
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.text = element_text(size = 17))
  
}

# DEFINE FUNCTION to generate fold change ranked list for GSEA -------------------
makeFCrankedList <- function(dataset_name, DE_res_file, outPath, analysis){
  RankedListType <- "FoldChange"
  
  ### Generate ranked list
  res <- readRDS(DE_res_file)
  res <- cbind(rownames(res), res)
  colnames(res) <- c("gene","baseMean","logFC","lfcSE","stat","p.value","pFDR")
  
  ### Generate ranked list
  colNameGenes <- "gene"
  RankedList <- res$logFC
  names(RankedList) <- res[,colNameGenes]
  if(length( which(is.na(RankedList)))>0){
    RankedList <- RankedList[-(which(is.na(RankedList)))]
  }
  RankedList <- sort(RankedList)
  RankedList <- as.data.frame(cbind(names(RankedList),RankedList))
  colnames(RankedList) <- c("GeneName","RankedList")
  write.table(file=paste0(outPath, "FoldChangeRanked_", analysis, "_",dataset_name,".rnk"), RankedList, col.names = T, row.names = T, sep="\t")
}

