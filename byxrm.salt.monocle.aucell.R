#segregants 
#names.df = read.delim("/Users/noahalexander/sc.chemostat.out/2rpm/filtered_feature_bc_matrix/features.tsv", header = F)

cds.t0 <- load_mm_data(mat_path = "/Users/noahalexander/seg_cellranger_report/t0/filtered_feature_bc_matrix/matrix.mtx", 
  feature_anno_path = "/Users/noahalexander/seg_cellranger_report/t0/filtered_feature_bc_matrix/features.tsv", 
  cell_anno_path = "/Users/noahalexander/seg_cellranger_report/t0/filtered_feature_bc_matrix/barcodes.tsv")


cds.t30 <- load_mm_data(mat_path = "/Users/noahalexander/seg_cellranger_report/t30/filtered_feature_bc_matrix/matrix.mtx", 
  feature_anno_path = "/Users/noahalexander/seg_cellranger_report/t30/filtered_feature_bc_matrix/features.tsv", 
  cell_anno_path = "/Users/noahalexander/seg_cellranger_report/t30/filtered_feature_bc_matrix/barcodes.tsv")



cds.t0@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2
cds.t30@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2


#use this to generate cc assignments 
cds.t0 = pca.umap.cluster.cc(cds.t0)
cds.t30 = pca.umap.cluster.cc(cds.t30)

cds.t0 = mono.aucell(cds.t0)
cds.t30 = mono.aucell(cds.t30)

#resolution to generate 8 clusters/to have a bit more than the existing set of discrete cc labels 
cds.t30 = cluster_cells(cds.t30, resolution = 0.0009)
cds.t30@colData$clusters = cds.t30@clusters$UMAP$clusters
plot_cells(cds.t30,  label_cell_groups = F, cell_size = 1)

#see slides for context regareding assignment of monocle clusters to cc phase discrete classifications 
#here I am just assigning alpha as one of the cc phases for convenience but obviously it is not a cc phase

#cln2 and htb2 expression is very close in t30 so i am subsetting to a few relevant clusters to 
#hopefully redo embedding with only these cells and get something a little clearer..
#cds.t30.sub = cds.t30[,cds.t30@colData$clusters %in% c(1,8,2)]
#cds.t30.sub = pca.umap.cluster.cc(cds.t30.sub)
#cds.t30.sub = cluster_cells(cds.t30.sub, resolution = 0.0005)


cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters ==1] = "S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==2] = "S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==3] = "M/G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==4] = "G1/S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==5] = "M/G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==6] = "M/G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==7] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==8] = "S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==9] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==10] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==11] = "G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==12] = "M/G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==13] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==14] = "MATALPHA"




cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters ==1] = "G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==2] = "S"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==3] = "M/G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==4] = "G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==5] = "M/G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==6] = "M/G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==7] = "G2/M"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==8] = "G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==9] = "G2/M"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==10] = "G2/M"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==11] = "G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==12] = "M/G1"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==13] = "G2/M"
cds.t30@colData$crude_cc_phase[cds.t30@clusters$UMAP$clusters  ==14] = "MATALPHA"


#used default clustering resolution for t0 
cds.t0@colData$clusters = cds.t0@clusters$UMAP$clusters

#filter out alpha cluster 
#cds.t0 = cds.t0[,cds.t0@colData$cluster %in% c(1,2,3,4,5,6)]
#cds.t30 = cds.t30[,cds.t30@colData$cluster %in% c(1,2,3,4)]

#see slides for context regareding assignment of monocle clusters to cc phase discrete classifications 
#here I am just assigning alpha as one of the cc phases for convenience but obviously it is not a cc phase

cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters ==1] = "G2/M"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==2] = "G2/M"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==3] = "G1/S"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==4] = "M/G1"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==5] = "S"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==6] = "G1"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==7] = "MATALPHA"

#nothing getting compressed for t0
cds.t0@colData$crude_cc_phase = cds.t0@colData$cc_phase




#visualize cc and esr markers to hopefully 

cds.t0 = pca.umap.cluster.cc(cds.t0)
cds.t30 = pca.umap.cluster.cc(cds.t30)




#########joint 



cds = combine_cds(list(cds.t0, cds.t30))
cds@colData$sample = as.factor(cds@colData$sample)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2


cds = pca.umap.cluster.cc(cds)

cds= pca.umap.cluster.esr(cds)

plot_cells(cds)

cds@colData$named_dataset[cds@colData$sample ==1] = "NaCl_0.7M_t0_3051_rep1"
cds@colData$named_dataset[cds@colData$sample ==2] = "NaCl_0.7M_t30_3051_rep1"


#making df for josh 
df = as.data.frame(cds@colData)
df$cluster = NULL
df$sample = NULL
df$cell_name = df$cell
df$cell = NULL

write.table(df, sep = "\t", quote = F, file = "NaCl_0.7M_3051_rep1.annotations.tsv", row.names = F, col.names = colnames(df))
