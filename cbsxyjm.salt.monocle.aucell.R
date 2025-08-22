names.df = read.delim("/Users/noahalexander/3004.salt/t0/filtered_feature_bc_matrix/features.tsv", header = F)


cds.t0 <- load_mm_data(mat_path = "/Users/noahalexander/3004.salt/t0/filtered_feature_bc_matrix/matrix.mtx.gz", 
                       feature_anno_path = "/Users/noahalexander/3004.salt/t0/filtered_feature_bc_matrix/features.tsv", 
                       cell_anno_path = "/Users/noahalexander/3004.salt/t0/filtered_feature_bc_matrix/barcodes.tsv")


cds.t30 <- load_mm_data(mat_path = "/Users/noahalexander/3004.salt/t30/filtered_feature_bc_matrix/matrix.mtx.gz", 
                        feature_anno_path = "/Users/noahalexander/3004.salt/t30/filtered_feature_bc_matrix/features.tsv.gz", 
                        cell_anno_path = "/Users/noahalexander/3004.salt/t30/filtered_feature_bc_matrix/barcodes.tsv")





cds.t0@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2
cds.t30@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2


#use this to generate cc assignments 
cds.t0 = pca.umap.cluster.cc(cds.t0)
cds.t30 = pca.umap.cluster.cc(cds.t30)

cds.t0 = mono.aucell(cds.t0)
cds.t30 = mono.aucell(cds.t30)


##
cds.t0 = cluster_cells(cds.t0, resolution = 0.00015)



#resolution to generate 8 clusters/to have a bit more than the existing set of discrete cc labels 
cds.t30 = cluster_cells(cds.t30, resolution = 0.00015)
cds.t30@colData$clusters = cds.t30@clusters$UMAP$clusters
plot_cells(cds.t30,  label_cell_groups = F )




cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters ==1] = "S"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==2] = "G2/M"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==3] = "G1/S"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==4] = "UNCLEAR"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==5] = "M/G1"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==6] = "G2/M"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==7] = "G2/M"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==8] = "G1"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==9] = "G2/M"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==10] = "UNCLEAR"
cds.t0@colData$cc_phase[cds.t0@clusters$UMAP$clusters  ==11] = "MATALPHA"



cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters ==1] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==2] = "G1/S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==3] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==4] = "UNCLEAR"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==5] = "S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==6] = "UNCLEAR"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==7] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==8] = "M/G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==9] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==10] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==11] = "S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==12] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==13] = "S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==14] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==15] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==16] = "G2/M"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==17] = "M/G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==18] = "S"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==19] = "G1"
cds.t30@colData$cc_phase[cds.t30@clusters$UMAP$clusters  ==20] = "MATALPHA"






cds = combine_cds(list(cds.t0, cds.t30))
cds@colData$sample = as.factor(cds@colData$sample)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2


cds = pca.umap.cluster.cc(cds)

#cds= pca.umap.cluster.esr(cds)

#plot_cells(cds)
cds = mono.aucell(cds)

cds@colData$named_dataset[cds@colData$sample ==1] = "NaCl_0.7M_3004_rep1_t0"
cds@colData$named_dataset[cds@colData$sample ==2] = "NaCl_0.7M_3004_rep1_t30"


#making df for josh 
df = as.data.frame(cds@colData)
df$seurat_cluster = cds@clusters$UMAP$clusters
df$cell_cycle = df$cc_phase
df$cell_name = df$cell
df$clusters = NULL
df$cc_phase = NULL
df$cluster = NULL
df$sample = NULL
df$cell = NULL

write.table(df, sep = "\t", quote = F, file = "NaCl_0.7M_3004_rep1.annotations.20230910.tsv", row.names = F, col.names = colnames(df))
