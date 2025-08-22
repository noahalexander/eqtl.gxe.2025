library(dplyr)
library(AUCell)

names.df = read.delim("/Users/noahalexander/3004.salt/t0/filtered_feature_bc_matrix/features.tsv", header = F)
#names.df = read.delim("/Users/noahalexander/sc.chemostat.out/2rpm/filtered_feature_bc_matrix/features.tsv", header = F)

cc_markers = c("SPO12", "SWI5",  "CLB1" , "FKH1" ,"HTB2", "CLN2" , "CTS1" , "DSE4",  "DSE2" , "PIR1"  ,"ASH1" )



pca.umap.cluster <- function(cds.input){
  cds = preprocess_cds(cds.input)
  cds = reduce_dimension(cds)
  cds = cluster_cells(cds)
  cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2
  cds@colData$log2numi = log(cds@colData$n.umi, base = 2)

  return(cds)
}


###same but with cc focus for pca
pca.umap.cluster.cc <- function(cds.input){
  
  cds = preprocess_cds(cds.input, use_genes = intersect(rownames(cds.input), botstein.naming$ORF))
  cds = reduce_dimension(cds)
  cds = cluster_cells(cds)
  cds@colData$log2numi = log(cds@colData$n.umi, base = 2)
  cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2
  cds@colData$clusters = cds@clusters$UMAP$clusters
  
return(cds)
}


pca.umap.cluster.esr <- function(cds.input){
  df = read.delim("pbio.2004050.csv", sep = ",")
  cds = preprocess_cds(cds.input, use_genes = intersect(rownames(cds.input), df$UID))
  cds = reduce_dimension(cds)
  cds@colData$clusters = cds@clusters$UMAP$clusters
  cds@colData$log2numi = log(cds@colData$n.umi, base = 2)
  cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2
  return(cds)
}


mono.markers <- function(cds.input){
marker_test_res <- top_markers(cds.input, group_cells_by="cluster", cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.05) %>%
  group_by(cell_group) %>%
  top_n(10, specificity)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
plot_genes_by_group(cds.input,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)

return(top_specific_marker_ids)
}



mono.aucell <- function(cds.input) {
  df = read.delim("pbio.2004050.csv", sep = ",")
  unique(df$Gene.Group)
  df.rp = subset(df, df$Gene.Group == "RP cluster")
  df.iesr = subset(df, df$Gene.Group == "iESR cluster")
  df.ribi = subset(df, df$Gene.Group == "RiBi (originally called PAC) cluster")
  #hot1 = read.csv("Downloads/groupbytf_1415887615.csv", sep = "\t")
  
  cds = cds.input
  #rownames(cds) = names.df$V1
  geneSets <- list(rp_genes=intersect(rownames( cds) , df.rp$UID),
    iesr_genes=intersect(rownames( cds) ,df.iesr$UID),
    ribi_genes=intersect(rownames(cds), df.ribi$UID))
  
  exprMatrix = normalized_counts(cds)
  cells_rankings <- AUCell_buildRankings(exprMatrix)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=2, assign=TRUE)
  
  
  df = as.data.frame(cells_AUC@assays@data$AUC)
  df = as.data.frame(t(df))
  
  cds@colData$rp.aucell = df$rp_genes
  cds@colData$iesr.aucell = df$iesr_genes
  cds@colData$ribi.aucell = df$ribi_genes
  
  return(cds)
}



mono.aucell.cc <- function(cds.input) {
  df = read.csv("botstein.cc.20240327.csv")
  dim(df)
  df = subset(df, df$Aggregate.Score > 5)
  dim(df)
  unique(df$Peak)
  df.g2m = subset(df, df$Peak == "G2/M")
  df.g1 = subset(df, df$Peak == "G1")
  df.s = subset(df, df$Peak == "S")
  df.sg2 = subset(df, df$Peak == "S/G2")
  df.mg1 = subset(df, df$Peak == "M/G1")
  
  #hot1 = read.csv("Downloads/groupbytf_1415887615.csv", sep = "\t")
  cds = cds.input
  #rownames(cds) = names.df$V1
  geneSets <- list(g2m.genes=intersect(rownames( cds) , df.g2m$ORF),
    g1.genes=intersect(rownames( cds) ,df.g1$ORF),
    s.genes=intersect(rownames(cds), df.s$ORF),
    sg2.genes=intersect(rownames(cds), df.sg2$ORF),
    mg1.genes= intersect(rownames(cds), df.mg1$ORF))
  
  exprMatrix = normalized_counts(cds)
  cells_rankings <- AUCell_buildRankings(exprMatrix)
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, nCores=2, assign=TRUE)
  
  
  df = as.data.frame(cells_AUC@assays@data$AUC)
  df = as.data.frame(t(df))
  
  cds@colData$g2m.aucell = df$g2m.genes
  cds@colData$g1.aucell = df$g1.genes
  cds@colData$s.aucell = df$s.genes
  cds@colData$sg2.aucell = df$sg2.genes
  cds@colData$mg1.aucell = df$mg1.genes
  
  return(cds)
}
