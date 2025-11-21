# MSigDB hallmarks enrichment with fgsea -------------------------------------------------------
## With signed p-values empty output, with log2FC too many ties
GSEA_MSigDBHallmark <- function(bulkType, analysis, res, outPath){
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
  library(ggplot2)
  options(digits=15)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
  msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
  msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)

  FoldChange_ranked <- res$log2FoldChange
  names(FoldChange_ranked) <- rownames(as.data.frame(res))
  if(length( which(is.na(FoldChange_ranked)))>0){
    FoldChange_ranked <- FoldChange_ranked[-(which(is.na(FoldChange_ranked)))]
  }
  write.table(file=paste0(outPath, "FoldChangeRanked_", analysis, "_",bulkType,".rnk"), FoldChange_ranked, col.names = F, row.names = T, sep="\t")

  fgseaRes_FC_H <- fgsea(pathways = msigdbr_list_H_nr, 
                       stats    = FoldChange_ranked)
  #saveRDS(file=paste0(outPath, "fgsea_msigdbHallmark_FoldChangeRanked_", analysis, "_",bulkType,".RDS"), fgseaRes_FC_H)
  fgseaRes_FC_H$padj <- as.numeric(as.vector(fgseaRes_FC_H$padj))
  
  for(pval_cutoff in c(0.05, 0.1)){
    png(paste0(outPath, "fgsea_msigdbHallmark_FC_", analysis, "_",bulkType,"_RawP",pval_cutoff,".png"), width=1000, height = 400)
    p1 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),], aes(x=NES, 
                                                                         y=pathway, 
                                                                         colour=pval, 
                                                                         size=size)) +
    geom_point() + geom_vline(xintercept=0) +
    expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
    labs(x="Normalized Enrichment Score", y="pathway", colour="p value", size="N genes") +
    ggtitle(paste0("enriched pathways raw p-value < ",pval_cutoff))
    p2 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),],aes(x=ES, 
                                                                        y=pathway, 
                                                                        colour=padj, 
                                                                        size=size)) +
    geom_point() + geom_vline(xintercept=0) +
    expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
    labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") 
    print(grid.arrange(p1,p2, ncol=2))
    dev.off()
    plotOut <- grid.arrange(p1,p2, ncol=2)
  }
  return(fgseaRes_FC_H)
  write.table(file=paste0(outPath, "fgsea_msigdbHallmark_FC_", analysis, "_",bulkType,"_RawP005.txt"), as.data.frame(fgseaRes_FC_H[fgseaRes_FC_H$pval <0.05,c(1:3,5:7)]), col.names = T, row.names = F, quote=F)
}


# MSigDB GO C5 enrichment with fgsea -------------------------------------------------------
## With signed p-values empty output, with log2FC too many ties
GSEA_MSigDBGeneOntologyC5 <- function(bulkType, analysis, res, outPath){
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
  library(ggplot2)
  options(digits=15)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "C5")
  msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
  msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)
  
  FoldChange_ranked <- res$log2FoldChange
  names(FoldChange_ranked) <- rownames(as.data.frame(res))
  if(length( which(is.na(FoldChange_ranked)))>0){
    FoldChange_ranked <- FoldChange_ranked[-(which(is.na(FoldChange_ranked)))]
  }
  write.table(file=paste0(outPath, "FoldChangeRanked_", analysis, "_",bulkType,".rnk"), FoldChange_ranked, col.names = F, row.names = T, sep="\t")
  
  fgseaRes_FC_H <- fgsea(pathways = msigdbr_list_H_nr, 
                         stats    = FoldChange_ranked)
  #saveRDS(file=paste0(outPath, "fgsea_msigdbHallmark_FoldChangeRanked_", analysis, "_",bulkType,".RDS"), fgseaRes_FC_H)
  fgseaRes_FC_H$padj <- as.numeric(as.vector(fgseaRes_FC_H$padj))
  
  for(pval_cutoff in c(0.05, 0.1)){
    png(paste0(outPath, "fgsea_msigdbGOC5_FC_", analysis, "_",bulkType,"_RawP",pval_cutoff,".png"), width=1000, height = 400)
    p1 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),], aes(x=NES, 
                                                                                    y=pathway, 
                                                                                    colour=pval, 
                                                                                    size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Normalized Enrichment Score", y="pathway", colour="p value", size="N genes") +
      ggtitle(paste0("enriched pathways raw p-value < ",pval_cutoff))
    p2 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),],aes(x=ES, 
                                                                                   y=pathway, 
                                                                                   colour=padj, 
                                                                                   size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") 
    print(grid.arrange(p1,p2, ncol=2))
    dev.off()
    plotOut <- grid.arrange(p1,p2, ncol=2)
  }
  
  padj_cutoff <- 0.1
  png(paste0(outPath, "fgsea_msigdbGOC5_FC_", analysis, "_",bulkType,"_padj",padj_cutoff,".png"), width=1000, height = 400)
  print(ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$padj<padj_cutoff)),], aes(x=NES, 
                                                                            y=pathway, 
                                                                            colour=padj, 
                                                                            size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Normalized Enrichment Score", y="pathway", colour="padj", size="N genes") +
      ggtitle(paste0("enriched pathways adj p-value < ",padj_cutoff)))
  dev.off()
  
  top15  <- head(fgseaRes_FC_H[order(fgseaRes_FC_H$padj),], n=15)
  png(paste0(outPath, "fgsea_msigdbGOC5_FC_", analysis, "_",bulkType,"_Top15adjPval.png"), width=1000, height = 400)
  print(ggplot(data=top15, aes(x=ES, y=pathway, colour=padj, size=size)) +
    geom_point() + geom_vline(xintercept=0) +
    expand_limits(x=(min(top15$NES))) + theme_bw() + 
    ggtitle(paste0(bulkType, " ", analysis, " top 15 GO MSigC5")) +
    labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") )
  dev.off()
  for(type in c("GOBP","GOCC","GOMF")){
    res <- fgseaRes_FC_H[grep(type,fgseaRes_FC_H$pathway),]
    top15  <- head(res[order(res$padj),], n=15)
    png(paste0(outPath, "fgsea_msigdbGOC5", type, "_FC_", analysis, "_",bulkType,"_Top15adjPval.png"), width=1000, height = 400)
    print(ggplot(data=top15, aes(x=ES, y=pathway, colour=padj, size=size)) +
          geom_point() + geom_vline(xintercept=0) +
          expand_limits(x=(min(top15$NES))) + theme_bw() + 
          ggtitle(paste0(bulkType, " ", analysis, " top 15 GO MSigC5")) +
          labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") )
    dev.off()
  }
  return(fgseaRes_FC_H)
  write.table(file=paste0(outPath, "fgsea_msigdbGOC5_FC_", analysis, "_",bulkType,"_RawP005.txt"), as.data.frame(fgseaRes_FC_H[fgseaRes_FC_H$pval <0.05,c(1:3,5:7)]), col.names = T, row.names = F, quote=F)
}

# VERSION 2: MSigDB hallmarks enrichment with fgsea -------------------------------------------------------
## Changes of this second version: i. the possibility to set type of gene identifier, ii. input is ranked list directly
## With signed p-values empty output, with log2FC too many ties
GSEA_MSigDBHallmark_v2 <- function(bulkType, analysis, RankedList, RankedListType, outPath, geneIdType){
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
  library(ggplot2)
  options(digits=15)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")
  if(geneIdType=="gene_symbol"){msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)}
  if(geneIdType=="ensembl_id"){msigdbr_list_H = split(x = h_gene_sets$ensembl_gene, f = h_gene_sets$gs_name)}
  
  msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)
  
  fgseaRes_FC_H <- fgsea(pathways = msigdbr_list_H_nr, 
                         stats    = RankedList)
  
  fgseaRes_FC_H$padj <- as.numeric(as.vector(fgseaRes_FC_H$padj))
  
  for(pval_cutoff in c(0.05, 0.1)){
    png(paste0(outPath, "fgsea_msigdbHallmark_", analysis, "_",bulkType,"_",RankedListType,"_RawP",pval_cutoff,".png"), width=1000, height = 400)
    p1 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),], aes(x=NES, 
                                                                                    y=pathway, 
                                                                                    #colour=pval, 
                                                                                    size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Normalized Enrichment Score", y="pathway", colour="p value", size="N genes") +
      ggtitle(paste0("enriched pathways raw p-value < ",pval_cutoff))
    p2 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),],aes(x=ES, 
                                                                                   y=pathway, 
                                                                                   colour=padj, 
                                                                                   size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") 
    print(grid.arrange(p1,p2, ncol=2))
    dev.off()
    plotOut <- grid.arrange(p1,p2, ncol=2)
  }
  
  padj_cutoff <- 0.1
  png(paste0(outPath, "fgsea_msigdbHallmark_", analysis, "_",bulkType,"_",RankedListType,"_padj",padj_cutoff,".png"), width=1000, height = 400)
  print(ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$padj<padj_cutoff)),], aes(x=NES, 
                                                                                  y=pathway, 
                                                                                  colour=padj, 
                                                                                  size=size)) +
    geom_point() + geom_vline(xintercept=0) +
    expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
    labs(x="Normalized Enrichment Score", y="pathway", colour="adj pvalue", size="N genes") +
    ggtitle(paste0("enriched pathways adj p-value < ",padj_cutoff)))
  dev.off()
  
  return(fgseaRes_FC_H)
  write.table(file=paste0(outPath, "fgsea_msigdbHallmark_", analysis, "_",bulkType,"_",RankedListType,"_RawP005.txt"), as.data.frame(fgseaRes_FC_H[fgseaRes_FC_H$pval <0.05,c(1:3,5:7)]), col.names = T, row.names = F, quote=F)
}


# VERSION2 MSigDB GO C5 enrichment with fgsea -------------------------------------------------------
## The only change of this second version is the possibility to set type of gene identifier
## With signed p-values empty output, with log2FC too many ties
GSEA_MSigDBGeneOntologyC5_v2 <- function(bulkType, analysis, RankedList, RankedListType, outPath, geneIdType){
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
  library(ggplot2)
  options(digits=15)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "C5")
  if(geneIdType=="gene_symbol"){msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)}
  if(geneIdType=="ensembl_id"){msigdbr_list_H = split(x = h_gene_sets$ensembl_gene, f = h_gene_sets$gs_name)}
  msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)
  
  fgseaRes_FC_H <- fgsea(pathways = msigdbr_list_H_nr, 
                         stats    = RankedList)
  fgseaRes_FC_H$padj <- as.numeric(as.vector(fgseaRes_FC_H$padj))
  
  for(pval_cutoff in c(0.05, 0.1)){
    png(paste0(outPath, "fgsea_msigdbGOC5_", analysis, "_",bulkType,"_",RankedListType,"_RawP",pval_cutoff,".png"), width=1000, height = 400)
    p1 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),], aes(x=NES, 
                                                                                    y=pathway, 
                                                                                    colour=pval, 
                                                                                    size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Normalized Enrichment Score", y="pathway", colour="p value", size="N genes") +
      ggtitle(paste0("enriched pathways raw p-value < ",pval_cutoff))
    p2 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),],aes(x=ES, 
                                                                                   y=pathway, 
                                                                                   colour=padj, 
                                                                                   size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") 
    print(grid.arrange(p1,p2, ncol=2))
    dev.off()
    plotOut <- grid.arrange(p1,p2, ncol=2)
  }
  
  padj_cutoff <- 0.1
  png(paste0(outPath, "fgsea_msigdbGOC5_", analysis, "_",bulkType,"_",RankedListType,"_padj",padj_cutoff,".png"), width=1000, height = 400)
  print(ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$padj<padj_cutoff)),], aes(x=NES, 
                                                                                  y=pathway, 
                                                                                  colour=padj, 
                                                                                  size=size)) +
    geom_point() + geom_vline(xintercept=0) +
    expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
    labs(x="Normalized Enrichment Score", y="pathway", colour="padj", size="N genes") +
    ggtitle(paste0("enriched pathways adj p-value < ",padj_cutoff)))
  dev.off()
  
  #top15  <- head(fgseaRes_FC_H[order(fgseaRes_FC_H$padj),], n=15)
  top15  <- head(fgseaRes_FC_H[order(fgseaRes_FC_H$padj),], n=40)
  #png(paste0(outPath, "fgsea_msigdbGOC5_", analysis, "_",bulkType,"_",RankedListType,"_Top15adjPval.png"), width=1000, height = 400)
  png(paste0(outPath, "fgsea_msigdbGOC5_", analysis, "_",bulkType,"_",RankedListType,"_Top40adjPval.png"), width=1000, height = 400)
  print(ggplot(data=top15, aes(x=ES, y=pathway, colour=padj, size=size)) +
          geom_point() + geom_vline(xintercept=0) +
          expand_limits(x=(min(top15$NES))) + theme_bw() + 
          ggtitle(paste0(bulkType, " ", analysis, " top 15 GO MSigC5")) +
          #ggtitle(paste0(bulkType, " ", analysis, " top 40 GO MSigC5")) +
          labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") )
  dev.off()
  for(type in c("GOBP","GOCC","GOMF")){
    res <- fgseaRes_FC_H[grep(type,fgseaRes_FC_H$pathway),]
    #top15  <- head(res[order(res$padj),], n=15)
    top40  <- head(res[order(res$padj),], n=40)
    #png(paste0(outPath, "fgsea_msigdbGOC5", type, "_FC_", analysis, "_",bulkType,"_",RankedListType,"_Top15adjPval.png"), width=1000, height = 400)
    png(paste0(outPath, "fgsea_msigdbGOC5", type, "_FC_", analysis, "_",bulkType,"_",RankedListType,"_Top40adjPval.png"), width=1000, height = 400)
    print(ggplot(data=top15, aes(x=ES, y=pathway, colour=padj, size=size)) +
            geom_point() + geom_vline(xintercept=0) +
            expand_limits(x=(min(top15$NES))) + theme_bw() + 
            ggtitle(paste0(bulkType, " ", analysis, " top 15 GO MSigC5")) +
            ggtitle(paste0(bulkType, " ", analysis, " top 40 GO MSigC5")) +
            labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") )
    dev.off()
  }
  return(fgseaRes_FC_H)
  write.table(file=paste0(outPath, "fgsea_msigdbGOC5_", analysis, "_",bulkType,"_",RankedListType,"_RawP005.txt"), as.data.frame(fgseaRes_FC_H[fgseaRes_FC_H$pval <0.05,c(1:3,5:7)]), col.names = T, row.names = F, quote=F)
}


# MSigDB Curated Gene Sets, Canonical Pathways (C2:CP) enrichment with fgsea -------------------------------------------------------
GSEA_MSigDBGeneOntologyC2CP <- function(bulkType, analysis, RankedList, RankedListType, outPath, geneIdType){
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
  library(ggplot2)
  options(digits=15)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory="CP")
  if(geneIdType=="gene_symbol"){msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)}
  if(geneIdType=="ensembl_id"){msigdbr_list_H = split(x = h_gene_sets$ensembl_gene, f = h_gene_sets$gs_name)}
  msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)
  
  
  fgseaRes_FC_H <- fgsea(pathways = msigdbr_list_H_nr, 
                         stats    = RankedList)
  #saveRDS(file=paste0(outPath, "fgsea_msigdbHallmark_FoldChangeRanked_", analysis, "_",bulkType,".RDS"), fgseaRes_FC_H)
  fgseaRes_FC_H$padj <- as.numeric(as.vector(fgseaRes_FC_H$padj))
  
  for(pval_cutoff in c(0.05, 0.1)){
    png(paste0(outPath, "fgsea_msigdbC2CP_", analysis, "_",bulkType,"_",RankedListType,"_RawP",pval_cutoff,".png"), width=1000, height = 400)
    p1 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),], aes(x=NES, 
                                                                                    y=pathway, 
                                                                                    colour=pval, 
                                                                                    size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Normalized Enrichment Score", y="pathway", colour="p value", size="N genes") +
      ggtitle(paste0("enriched pathways raw p-value < ",pval_cutoff))
    p2 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),],aes(x=ES, 
                                                                                   y=pathway, 
                                                                                   colour=padj, 
                                                                                   size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") 
    print(grid.arrange(p1,p2, ncol=2))
    dev.off()
    plotOut <- grid.arrange(p1,p2, ncol=2)
  }
  top15  <- head(fgseaRes_FC_H[order(fgseaRes_FC_H$padj),], n=15)
  png(paste0(outPath, "fgsea_msigdbC2CP_", analysis, "_",bulkType,"_",RankedListType,"_Top15adjPval.png"), width=1000, height = 400)
  print(ggplot(data=top15, aes(x=ES, y=pathway, colour=padj, size=size)) +
          geom_point() + geom_vline(xintercept=0) +
          expand_limits(x=(min(top15$NES))) + theme_bw() + 
          ggtitle(paste0(bulkType, " ", analysis, " top 15 GO MSigC2CP")) +
          labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") )
  dev.off()
  
  return(fgseaRes_FC_H)
  write.table(file=paste0(outPath, "fgsea_msigdbC2CP_", analysis, "_",bulkType,"_",RankedListType,"_RawP005.txt"), as.data.frame(fgseaRes_FC_H[fgseaRes_FC_H$pval <0.05,c(1:3,5:7)]), col.names = T, row.names = F, quote=F)
}


# GSEA on Transcription Factor Binding Sites (msigdb C3) ------------------
GSEA_MSigDB_TFBS <- function(bulkType, analysis, RankedList, RankedListType, outPath, geneIdType){
  library(msigdbr)
  library(fgsea)
  library(gridExtra)
  library(ggplot2)
  options(digits=15)
  h_gene_sets = msigdbr(species = "Homo sapiens", category = "C3")
  msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
  if(geneIdType=="gene_symbol"){msigdbr_list_H = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)}
  if(geneIdType=="ensembl_id"){msigdbr_list_H = split(x = h_gene_sets$ensembl_gene, f = h_gene_sets$gs_name)}
  
  msigdbr_list_H_nr <- lapply(msigdbr_list_H, unique)
  
  fgseaRes_FC_H <- fgsea(pathways = msigdbr_list_H_nr, 
                         stats    = RankedList)
  
  fgseaRes_FC_H$padj <- as.numeric(as.vector(fgseaRes_FC_H$padj))
  
  for(pval_cutoff in c(0.05, 0.1)){
    png(paste0(outPath, "fgsea_msigdbC3TFBS_", analysis, "_",bulkType,"_",RankedListType,"_RawP",pval_cutoff,".png"), width=1000, height = 400)
    p1 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),], aes(x=NES, 
                                                                                    y=pathway, 
                                                                                    #colour=pval, 
                                                                                    size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Normalized Enrichment Score", y="pathway", colour="p value", size="N genes") +
      ggtitle(paste0("enriched pathways raw p-value < ",pval_cutoff))
    p2 <- ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$pval<pval_cutoff)),],aes(x=ES, 
                                                                                   y=pathway, 
                                                                                   colour=padj, 
                                                                                   size=size)) +
      geom_point() + geom_vline(xintercept=0) +
      expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
      labs(x="Enrichment Score", y="pathway", colour="adjusted p value", size="N genes") 
    print(grid.arrange(p1,p2, ncol=2))
    dev.off()
    plotOut <- grid.arrange(p1,p2, ncol=2)
  }
  
  padj_cutoff <- 0.1
  png(paste0(outPath, "fgsea_msigdbC3TFBS_", analysis, "_",bulkType,"_",RankedListType,"_padj",padj_cutoff,".png"), width=1000, height = 400)
  print(ggplot(data= fgseaRes_FC_H[which((fgseaRes_FC_H$padj<padj_cutoff)),], aes(x=NES, 
                                                                                  y=pathway, 
                                                                                  colour=padj, 
                                                                                  size=size)) +
          geom_point() + geom_vline(xintercept=0) +
          expand_limits(x=(min(fgseaRes_FC_H$NES))) + theme_bw() +
          labs(x="Normalized Enrichment Score", y="pathway", colour="adj pvalue", size="N genes") +
          ggtitle(paste0("enriched pathways adj p-value < ",padj_cutoff)))
  dev.off()
  
  return(fgseaRes_FC_H)
  write.table(file=paste0(outPath, "fgsea_msigdbC3TFBS", analysis, "_",bulkType,"_",RankedListType,"_RawP005.txt"), as.data.frame(fgseaRes_FC_H[fgseaRes_FC_H$pval <0.05,c(1:3,5:7)]), col.names = T, row.names = F, quote=F)
}
# Volcano plots -----------------------------------------------------------
volcanoPlot <- function(analysis, bulkType, res, betaCatTargets){
  library("ggrepel")
  library("ggplot2")
  df <- as.data.frame(cbind(rownames(res), res))
  colnames(df)[1] <- "gene" 
  df <- df[-(which(is.na(df$pvalue))),]

  #df[df$pvalue==0,]$pvalue <- .Machine$double.xmin # i replace p-values of 0 (which would be Inf in log transformation) with lowest machine value
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
  df$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 (FC=1.5) and pvalue < 0.05
  df$diffexpressed[df$log2FoldChange > 0.6 & df$pvalue < 0.05] <- "UP"
  df$diffexpressed[df$log2FoldChange < -0.6 & df$pvalue < 0.05] <- "DOWN"
  df$label <- ""

  ### Plot labeling B-cat targets with p_adj<0.1
  targetsToLabel <- df[((df$pvalue<0.05)&(df$gene %in% betaCatTargets)),]$gene
  df[df$gene %in% targetsToLabel,]$label <- df[df$gene %in% targetsToLabel,]$gene
  
  plotOut1 <- ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=label)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
    geom_point(size = 1, alpha=0.5) +
    geom_label_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("#00AFBB", "grey", "red"), # to set the colours of our variable
                       labels = c("Downregulated (L2FC<-0.6)", "Not significant (pval>0.05)", "Upregulated (L2FC>0.6)")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)</p><br /><br />
    ggtitle(paste0(analysis,"_",bulkType,"; labels on B-catenin targets pval<0.05 (no L2FC cutoff)"))
 
  ### Plot labeling genes with raw -log10(p-value) > 5, abs log2 FC >= 0.6
  df$label1 <- ""
  #df[((abs(df$log2FoldChange) >= 6) & (-log10(df$pvalue)>5)),]$label1 <- df[((abs(df$log2FoldChange) >= 6) & (-log10(df$pvalue)>5)),]$gene
  df[((-log10(df$pvalue)>4)&(df$log2FoldChange < -5)),]$label1 <- df[((-log10(df$pvalue)>4)&(df$log2FoldChange < -5)),]$gene
  df[((-log10(df$pvalue)>4)&(df$log2FoldChange > 5)),]$label1 <- df[((-log10(df$pvalue)>4)&(df$log2FoldChange > 5)),]$gene

  plotOut2 <- ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=label1)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values = c("#00AFBB", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)</p><br /><br />
  ggtitle(paste0(analysis,"_",bulkType,"; labels on genes pval<0.0001, |L2FC| > 5"))
  
  return(list(BcatTargetsPadj01=plotOut1, RawMinusLogPval5_log2FC05=plotOut2))
}  

# Volcano plots public ACC data-----------------------------------------------------------
# Highligthing sign DE B-cat targets or highlighting independently of significance, 9 selected B-cat targets from B-cat activation
# signature (targets sign up in at least 3 of the 4 public ACC data in CTNNB1 mut vs wt comparisons).
# Need in input df with at least colupmn "logFC", "gene", "pFDR
volcanoPlots_BcatTargets <- function(outPath, res, dataset, analysis){
  NusseHerbstPDE2A <- scan(file="C:/Users/LF260934/Documents/WORKING_FOLDER/Results/ACC_Cochin_snRNAseq/MetadataAndCo/BCateninTargets_Nusse_Herbst2014_PDE2A_nr.txt", what=character())
  betacat_targets <- c(NusseHerbstPDE2A, "AFF3")
  betacat_targets <- betacat_targets[order(betacat_targets)]
  res$pFDR <- as.numeric(as.vector(res$pFDR))
  res$logFC <- as.numeric(as.vector(res$logFC))
  
  ### Plot labeling B-cat targets with p_adj<0.1
  targetsToLabel <- res[((res$pFDR<0.1)&(res$gene %in% betacat_targets)),]$gene
  res$label <- ""
  res[res$gene %in% targetsToLabel,]$label <- res[res$gene %in% targetsToLabel,]$gene
  res$col_label <- ifelse(res$label=="", "gray", "black")
  res$col_label <- factor(res$col_label, levels=c("gray","black"))
  plot1 <- ggplot(data = res, aes(x = logFC, y = -log10(pFDR), col = col_label, label=label)) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.1), col = "black", linetype = 'dashed') +
    geom_point(size = 1, alpha=0.5) +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("gray","black")) +
    ggtitle(paste0(dataset, ", ", analysis, " (label sign DE BcatTargets)")) + theme_bw()
  
  # Show selected B-cat targets independently of significance
  res$label1 <- ""
  selectionOnPublicData9genes <- c("AXIN2","AFF3","LEF1","FSCN1", "MYC", "SP5","PDE2A","TBX3","TNFRSF19") # target genes found to be most consistently up in CTNNB1 mut vs wt in Assie,Heaton,TCGA,(Lefevre, only qualitative)
  targetsToLabel <- res[(res$gene %in% selectionOnPublicData9genes),]$gene
  res[res$gene %in% targetsToLabel,]$label1 <- res[res$gene %in% targetsToLabel,]$gene
  res$col_label <- ifelse(res$label1=="", "gray", "black")
  res$col_label <- factor(res$col_label, levels=c("gray","black"))
  plot2 <- ggplot(data = res, aes(x = logFC, y = -log10(pFDR), col = col_label, label=label1)) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.1), col = "black", linetype = 'dashed') +
    geom_point(size = 1, alpha=0.5) +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("gray","black")) + # to set the colours of our variable
    ggtitle(paste0(dataset, ", ", analysis, " (label 9 BcatTargets signature)")) + theme_bw()
  
  return(list(BcatSignificant=plot1, NineSignatureTargets=plot2))
}

# Volcano plots: highhlight a selected list of genes (in input of the function) -----------------------------------------------------------
volcanoPlots_highlightGenes <- function(res, dataset, analysis, genes){
  res$pFDR <- as.numeric(as.vector(res$pFDR))
  res$logFC <- as.numeric(as.vector(res$logFC))
  
  # Show selected genes independently of significance
  res$label1 <- ""
  targetsToLabel <- res[(res$gene %in% genes),]$gene
  res[res$gene %in% targetsToLabel,]$label1 <- res[res$gene %in% targetsToLabel,]$gene
  res$col_label <- ifelse(res$label1=="", "gray", "black")
  res$col_label <- factor(res$col_label, levels=c("gray","black"))
  plot2 <- ggplot(data = res, aes(x = logFC, y = -log10(pFDR), col = col_label, label=label1)) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "black", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.1), col = "black", linetype = 'dashed') +
    geom_point(size = 1, alpha=0.5) +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("gray","black")) + # to set the colours of our variable
    ggtitle(paste0(dataset, ", ", analysis)) + theme_bw()
  
  return(plot2)
}