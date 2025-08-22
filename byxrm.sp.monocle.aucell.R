#segregants 
#names.df = read.delim("/Users/noahalexander/sc.chemostat.out/2rpm/filtered_feature_bc_matrix/features.tsv", header = F)

cds.t0 <- load_mm_data(mat_path = "/Users/noahalexander/q.seg.2023/t0/filtered_feature_bc_matrix/matrix.mtx", 
                       feature_anno_path = "/Users/noahalexander/q.seg.2023/t0/filtered_feature_bc_matrix/features.tsv", 
                       cell_anno_path = "/Users/noahalexander/q.seg.2023/t0/filtered_feature_bc_matrix/barcodes.tsv")


cds.t10 <- load_mm_data(mat_path = "/Users/noahalexander/q.seg.2023/t10/filtered_feature_bc_matrix/matrix.mtx", 
                        feature_anno_path = "/Users/noahalexander/q.seg.2023/t10/filtered_feature_bc_matrix/features.tsv", 
                        cell_anno_path = "/Users/noahalexander/q.seg.2023/t10/filtered_feature_bc_matrix/barcodes.tsv")



cds.t0@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2

cds.t10@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2





#use this to generate cc assignments 
cds.t0 = pca.umap.cluster.cc(cds.t0)
cds.t10 = pca.umap.cluster.cc(cds.t10)

#cds.t10 = pca.umap.cluster.esr(cds.t10)


cds.t0 = mono.aucell(cds.t0)
cds.t10 = mono.aucell(cds.t10)



###########

#good resultion for t0 segs
cds.t0 = cluster_cells(cds.t0, resolution = 0.0003)
plot_cells(cds.t0, cell_size = 1.5, label_cell_groups = T, group_label_size = 15)

top_specific_marker_ids = mono.markers(cds.t0)
plot_cells(cds.t0, cell_size = 1.5, label_cell_groups = T, group_label_size = 15)
plot_cells(cds.t0, cell_size = 1.5, genes = "GRE1")

plot_genes_by_group(cds.t0,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)


#DSP = deep stationary phase
#ssp = shallow tationary phase
#ssp_cc = shallow stationary phase with cc markers at a low level 
#pseudo_transporter = expression of high-affinity glu transporters and pseudohyphal growth markers
#extreme = small subset of unusal cells I am tagging as such without assigning a clear label 

cds.t0@colData$state[cds.t0@clusters$UMAP$clusters ==1] = "DSP"
cds.t0@colData$state[cds.t0@clusters$UMAP$clusters  ==2] = "SSP_CC"
cds.t0@colData$state[cds.t0@clusters$UMAP$clusters  ==3] = "DSP"
cds.t0@colData$state[cds.t0@clusters$UMAP$clusters  ==4] = "SSP"
cds.t0@colData$state[cds.t0@clusters$UMAP$clusters  ==5] = "pseudo_transporter"
cds.t0@colData$state[cds.t0@clusters$UMAP$clusters  ==6] = "DSP"
cds.t0@colData$state[cds.t0@clusters$UMAP$clusters  ==7] = "EXTREME"




###############

cds.t10 = cluster_cells(cds.t10, resolution = 0.0003)
plot_cells(cds.t10, cell_size = 1.5, label_cell_groups = T, group_label_size = 15)

top_specific_marker_ids = mono.markers(cds.t10)
plot_cells(cds.t10, reduction_method = 'UMAP', genes = cc_markers, cell_size = 2)

plot_genes_by_group(cds.t10,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)





cds.t10 = readRDS("cds.t10.nocc.RDS")

cds.t10.nf = cds.t10 


cds.t10@colData$state[cds.t10@clusters$UMAP$clusters ==1] = "dessication_rec"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==2] = "UNCLEAR_A"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==3] = "STRESS_ribi_off"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==4] = "STRESS_ribi_off"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==5] = "dessication_rec"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==6] = "M/G1"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==7] = "STRESS_ribi_on"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==8] = "UNCLEAR_B"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==9] = "UNCLEAR_A"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==10] = "UNCLEAR_B"
cds.t10@colData$state[cds.t10@clusters$UMAP$clusters  ==11] = "UNCLEAR_A"

cds.t10@colData$cc_phase = cds.t10@colData$state


#resolution to generate 8 clusters/to have a bit more than the existing set of discrete cc labels 
#cds.t10 = cluster_cells(cds.t10, resolution = 0.00065)
#cds.t10@colData$clusters = cds.t10@clusters$UMAP$clusters
#plot_cells(cds.t10,  label_cell_groups = F, cell_size = 1)

#cds.t0 = cluster_cells(cds.t0, resolution = 0.00065)
#cds.t0@colData$clusters = cds.t0@clusters$UMAP$clusters
#plot_cells(cds.t0,  label_cell_groups = F, cell_size = 1)

 #cds.t0 = readRDS("cds.q.segs.cds.t0.20230901.RDS")
 #cds.t10 = readRDS("cds.q.segs.cds.t10.20230901.RDS")

cds = combine_cds(list(cds.t0, cds.t10))
cds@colData$sample = as.factor(cds@colData$sample)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] = names.df$V2


cds@colData$condition[cds@colData$sample ==1] = "SP 1wk"
cds@colData$condition[cds@colData$sample ==2] = "SP 1wk + YPD 10min"



cds = preprocess_cds(cds)
cds = reduce_dimension(cds)

#subset alphas using counts vs umap since there isn't a clear cluster and mfalpha1/YPL187W expression is speckled 
df = normalized_counts(cds)
df= as.matrix(df)
df = as.data.frame(t(df))
df$cell_name = gsub('.{0,2}$', '', rownames(df))

dfa = subset(df, df$YPL187W > 0)


#subset cds object for those cells 
cds.alpha = cds[, which(colnames(cds) %in% rownames(dfa))]
cds@colData$state

#assign mfalphas based on presence in the 'extreme' cluster or simple alpha expression (many but not 41 of the extreme cells are alphas so I am adding the rest of the extremes since 41 isn't enough to do anything with and they are prob alphas )
cds@colData$state[cds@colData$cell %in% dfa$cell_name] = "MFALPHA"
cds@colData$state[cds@colData$state %in% "EXTREME"] = "MFALPHA"


cds = pca.umap.cluster.cc(cds)

cds= pca.umap.cluster.esr(cds)


plot_cells(cds)
plot_cells(cds, color_cells_by = "sample")
plot_cells(cds, color_cells_by = "log2numi")
plot_cells(cds, color_cells_by = "iesr.aucell")
plot_cells(cds, color_cells_by = "ribi.aucell")
plot_cells(cds, color_cells_by = "rp.aucell")










top_specific_marker_ids = mono.markers(cds)














cds.q = cds
cds.salt = readRDS("cds.salt.segs.20230820.RDS")
cds.salt@colData$condition[cds.salt@colData$sample ==1] = "Mid-log YPD"
cds.salt@colData$condition[cds.salt@colData$sample ==2] = "Mid-log YPD + 0.7M NaCl 30min"

cds.all = combine_cds(list(cds.q, cds.salt))
cds.all = pca.umap.cluster.cc(cds.all)


plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)





#######where cds is t0 and t10


cds@colData$named_dataset[cds@colData$sample ==1] = "SP_3051_rep1_t0"
cds@colData$named_dataset[cds@colData$sample ==2] = "SP_3051_rep1_t10"


#making df for josh 
df = as.data.frame(cds@colData)
df$seurat_cluster = cds@clusters$UMAP$clusters
df$cell_cycle = df$state
df$cell_name = df$cell
df$clusters = NULL
df$cc_phase = NULL
df$cluster = NULL
df$sample = NULL
df$cell = NULL

head(df)

write.table(df, sep = "\t", quote = F, file = "SP_3051_rep1.annotations.20230910.tsv", row.names = F, col.names = colnames(df))
