
#3051 nacl t0

df = read.delim("/Users/noahalexander/Downloads/NaCl_0.7M_3051_rep1.extended.annotations.tsv")
d0 = subset(df, df$named_dataset =="NaCl_0.7M_t0_3051_rep1")
d30 = subset(df, df$named_dataset =="NaCl_0.7M_t30_3051_rep1")

#susbet to avoid retaining old cc labels/simplify things 
d0 = data.frame(cell_name= d0$cell_name, ribi.aucell = d0$ribi.aucell, rp.aucell=d0$rp.aucell, iesr.aucell=d0$iesr.aucell)
d30 = data.frame(cell_name= d30$cell_name, ribi.aucell = d30$ribi.aucell, rp.aucell=d30$rp.aucell, iesr.aucell=d30$iesr.aucell)



pbmc = readRDS("3051.naclt0.seurat.mating.cclabels.RDS")
pbmc@meta.data$cell_cycle = pbmc@active.ident
pbmc=seurat.aucell(pbmc)
md.3051.t0 = pbmc@meta.data
md.3051.t0$cell_name = rownames(md.3051.t0)


#finalized t30 object where alphas have been removed
pbmc = readRDS("3051.naclt30.seurat.mating.cclabels.RDS")
pbmc@meta.data$cell_cycle = pbmc@active.ident
pbmc=seurat.aucell(pbmc)
md.3051.t30 = pbmc@meta.data
md.3051.t30$cell_name = rownames(md.3051.t30)


#read in unfiltered data (from filtered bc matrix but not further filtered outside of cellranger)
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t0/filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
#dim(pbmc)
#6572 16160
unfiltered.3051.naclt0.md = pbmc@meta.data
unfiltered.3051.naclt0.md$cell_name = rownames(unfiltered.3051.naclt0.md)

#first merge the unfiltered 3051 nacl t0 metadata with the aucell annotations from the seurat df.
#merge the assignments such that all existing ribi assignments go to the unfiltered pbmc metadata and N/As are added for 67 cells in seurat not in monocle 
result_merge.t0 <- merge(unfiltered.3051.naclt0.md, d0, by = "cell_name", all.x = TRUE)

#next, merge the seurat-based cell cycle labels from the df where alphas have already been removed 
#the cell_cycle column should contain discrete assignment from pca ananlysis in seurat as well as NAs for cells removed (here, just alphas)  
result_merge.test2 <- merge(result_merge.t0, md.3051.t0, by = "cell_name", all.x = TRUE)

#add named_data 
result_merge.test2$named_data = rep("NaCl_0.7M_t0_3051_rep1", length(rownames(result_merge.test2)))
result_merge.test2$n.umi = result_merge.test2$nCount_RNA.x
result_merge.test2$log2numi =  log(result_merge.test2$nCount_RNA.y, base = 2)

result_merge.test2$cell_cycle = as.character(result_merge.test2$cell_cycle)
result_merge.test2$cell_cycle[is.na(result_merge.test2$cell_cycle)] <- "MATALPHA"

result_merge.test2$orig.ident.x <- NULL
result_merge.test2$nCount_RNA.x <- NULL
result_merge.test2$orig.ident.y <- NULL
result_merge.test2$nCount_RNA.y <- NULL
result_merge.test2$nFeature_RNA.y <- NULL
result_merge.test2$nCount_SCT <- NULL
result_merge.test2$nFeature_SCT <- NULL
result_merge.test2$nFeature_RNA.x  <- NULL


#cc labels are remade without spaaces or slashes
g2m.no.mating.indices = which(dd$cell_cycle == "G2_M w/o mating")
g2m.mating.indices = which(dd$cell_cycle == "G2_M + mating")

mg1.mating.indices = which(dd$cell_cycle == "M_G1 + mating")
mg1.no.mating.indices = which(dd$cell_cycle == "M/G1 w/o mating")


result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, g2m.no.mating.indices, "G2_M_no_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, g2m.mating.indices, "G2_M_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, mg1.mating.indices, "M_G1_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, mg1.no.mating.indices, "M_G1_no_mating")


#result_merge.test2, the nacl t0 table, is ready and we make the same thing for the 3051 nacl t30

#nacl t30


pbmc = readRDS("3051.naclt30.cclabels.20240528.RDS")
pbmc@meta.data$cell_cycle = pbmc@active.ident
md.3051.t30 = pbmc@meta.data
md.3051.t30$cell_name = rownames(md.3051.t30)


#read in unfiltered data (from filtered bc matrix but not further filtered outside of cellranger)
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t30/outs//filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
#dim(pbmc)
#6572 16160
unfiltered.3051.naclt30.md = pbmc@meta.data
unfiltered.3051.naclt30.md$cell_name = rownames(unfiltered.3051.naclt30.md)

#first merge the unfiltered 3051 nacl t0 metadata with the aucell annotations from the seurat df.
#merge the assignments such that all existing ribi assignments go to the unfiltered pbmc metadata and N/As are added for 67 cells in seurat not in monocle 
result_merge.t30 <- merge(unfiltered.3051.naclt30.md, d30, by = "cell_name", all.x = TRUE)

#next, merge the seurat-based cell cycle labels from the df where alphas have already been removed 
#the cell_cycle column should contain discrete assignment from pca ananlysis in seurat as well as NAs for cells removed (here, just alphas)  
result_merge.2.30 <- merge(result_merge.t30, md.3051.t30, by = "cell_name", all.x = TRUE)

#add named_data 
result_merge.2.30$named_data = rep("NaCl_0.7M_t30_3051_rep1", length(rownames(result_merge.2.30 )))
result_merge.2.30$n.umi = result_merge.2.30$nCount_RNA.x
result_merge.2.30$log2numi =  log(result_merge.2.30$nCount_RNA.y, base = 2)

result_merge.2.30$cell_cycle = as.character(result_merge.2.30$cell_cycle)
result_merge.2.30$cell_cycle[is.na(result_merge.2.30$cell_cycle)] <- "MATALPHA"

result_merge.2.30$orig.ident.x <- NULL
result_merge.2.30$nCount_RNA.x <- NULL
result_merge.2.30$orig.ident.y <- NULL
result_merge.2.30$nCount_RNA.y <- NULL
result_merge.2.30$nFeature_RNA.y <- NULL
result_merge.2.30$nCount_SCT <- NULL
result_merge.2.30$nFeature_SCT <- NULL
result_merge.2.30$nFeature_RNA.x  <- NULL



#write.table(result_merge.2.30, file= "3051.naclt30.statelabels.20240528.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)
#write.table(result_merge.test2, file= "3051.naclt0.statelabels.20240528.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)

dim(result_merge.test2)
dim(result_merge.2.30)

head(result_merge.test2)
head(result_merge.2.30)


test.out = rbind.fill(result_merge.test2, result_merge.2.30)
saveRDS(test.out, file = "3051.0.7MNaCl.states.2.RDS")
write.table(test.out, file= "3051.nacl.statelabels.20240529.2.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)



################--------------3051 sp experiments


#3051 sp t0

df = read.delim("/Users/noahalexander/Downloads/SP_3051_rep1.annotations.20230910.tsv")
d0 = subset(df, df$named_dataset =="SP_3051_rep1_t0")
d10 = subset(df, df$named_dataset =="SP_3051_rep1_t10")

#susbet to avoid retaining old cc labels/simplify things 
d0 = data.frame(cell_name= d0$cell_name, ribi.aucell = d0$ribi.aucell, rp.aucell=d0$rp.aucell, iesr.aucell=d0$iesr.aucell)
d10 = data.frame(cell_name= d10$cell_name, ribi.aucell = d10$ribi.aucell, rp.aucell=d10$rp.aucell, iesr.aucell=d10$iesr.aucell)



pbmc = readRDS("3051.sp.t0.fiveclusters.RDS")
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)

sp.t0.state = c("I", "II", "III", "IV", "V", "VI")
names(sp.t0.state) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, sp.t0.state)
pbmc@meta.data$cell_cycle = pbmc@active.ident
md.3051.t0 = pbmc@meta.data
md.3051.t0$cell_name = rownames(md.3051.t0)


#finalized t10 object where alphas have been removed
pbmc = readRDS("3051.spt10.lowresclustering.RDS")
sp.t10.state = c("I", "II", "III", "IV", "V", "VI","VII")
names(sp.t10.state) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, sp.t10.state)
pbmc@meta.data$cell_cycle = pbmc@active.ident
md.3051.t10 = pbmc@meta.data
md.3051.t10$cell_name = rownames(md.3051.t10)


#read in unfiltered data (from filtered bc matrix but not further filtered outside of cellranger)
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/seurat_in/")
pbmc <- CreateSeuratObject(pbmc_data)
#dim(pbmc)
#6572 16160
unfiltered.3051.spt0.md = pbmc@meta.data
unfiltered.3051.spt0.md$cell_name = rownames(unfiltered.3051.spt0.md)

#first merge the unfiltered 3051 nacl t0 metadata with the aucell annotations from the seurat df.
#merge the assignments such that all existing ribi assignments go to the unfiltered pbmc metadata and N/As are added for 67 cells in seurat not in monocle 
result_merge.t0 <- merge(unfiltered.3051.spt0.md, d0, by = "cell_name", all.x = TRUE)

#next, merge the seurat-based cell cycle labels from the df where alphas have already been removed 
#the cell_cycle column should contain discrete assignment from pca ananlysis in seurat as well as NAs for cells removed (here, just alphas)  
result_merge.test2 <- merge(result_merge.t0, md.3051.t0, by = "cell_name", all.x = TRUE)

#add named_data 
result_merge.test2$named_data = rep("SP_3051_rep1_t0", length(rownames(result_merge.test2)))
result_merge.test2$n.umi = result_merge.test2$nCount_RNA.x
result_merge.test2$log2numi =  log(result_merge.test2$nCount_RNA.y, base = 2)

result_merge.test2$cell_cycle = as.character(result_merge.test2$cell_cycle)
result_merge.test2$cell_cycle[is.na(result_merge.test2$cell_cycle)] <- "MATALPHA"

result_merge.test2$orig.ident.x <- NULL
result_merge.test2$nCount_RNA.x <- NULL
result_merge.test2$orig.ident.y <- NULL
result_merge.test2$nCount_RNA.y <- NULL
result_merge.test2$nFeature_RNA.y <- NULL
result_merge.test2$nCount_SCT <- NULL
result_merge.test2$nFeature_SCT <- NULL
result_merge.test2$nFeature_RNA.x  <- NULL

result_merge.test2.spt0 = result_merge.test2

#result_merge.test2, the sp t0 table, is ready and we make the same thing for the 3051 sp t10

#sp t10



#read in unfiltered data (from filtered bc matrix but not further filtered outside of cellranger)
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t10//filtered_feature_bc_matrix/seurat_in/")
pbmc <- CreateSeuratObject(pbmc_data)
#dim(pbmc)
#6572 16160
unfiltered.3051.spt10.md = pbmc@meta.data
unfiltered.3051.spt10.md$cell_name = rownames(unfiltered.3051.spt10.md)

#first merge the unfiltered 3051 nacl t0 metadata with the aucell annotations from the seurat df.
#merge the assignments such that all existing ribi assignments go to the unfiltered pbmc metadata and N/As are added for 67 cells in seurat not in monocle 
result_merge.t10 <- merge(unfiltered.3051.spt10.md, d10, by = "cell_name", all.x = TRUE)

#next, merge the seurat-based cell cycle labels from the df where alphas have already been removed 
#the cell_cycle column should contain discrete assignment from pca ananlysis in seurat as well as NAs for cells removed (here, just alphas)  
result_merge.2.10 <- merge(result_merge.t10, md.3051.t10, by = "cell_name", all.x = TRUE)

#add named_data 
result_merge.2.10$named_data = rep("SP_3051_rep1_t10", length(rownames(result_merge.2.10 )))
result_merge.2.10$n.umi = result_merge.2.10$nCount_RNA.x
result_merge.2.10$log2numi =  log(result_merge.2.10$nCount_RNA.y, base = 2)

result_merge.2.10$cell_cycle = as.character(result_merge.2.10$cell_cycle)
result_merge.2.10$cell_cycle[is.na(result_merge.2.10$cell_cycle)] <- "MATALPHA"

result_merge.2.10$orig.ident.x <- NULL
result_merge.2.10$nCount_RNA.x <- NULL
result_merge.2.10$orig.ident.y <- NULL
result_merge.2.10$nCount_RNA.y <- NULL
result_merge.2.10$nFeature_RNA.y <- NULL
result_merge.2.10$nCount_SCT <- NULL
result_merge.2.10$nFeature_SCT <- NULL
result_merge.2.10$nFeature_RNA.x  <- NULL



#write.table(result_merge.2.30, file= "3051.naclt30.statelabels.20240528.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)
#write.table(result_merge.test2, file= "3051.naclt0.statelabels.20240528.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)

dim(result_merge.test2)
dim(result_merge.2.10)

head(result_merge.test2)
head(result_merge.2.10)


test.out = rbind.fill(result_merge.test2, result_merge.2.10)

saveRDS(test.out, file = "3051.SP.states.2.RDS")
write.table(test.out, file= "3051.sp.statelabels.20240529.2.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)




#######-----3004 nacl 




#3004 nacl t0

#use this version, not the one without the (1)
df = read.delim("/Users/noahalexander/Downloads/NaCl_0.7M_3004_rep1.annotations (1).tsv")
d0 = subset(df, df$named_dataset =="NaCl_0.7M_3004_rep1_t0")
d30 = subset(df, df$named_dataset =="NaCl_0.7M_3004_rep1_t30")

#susbet to avoid retaining old cc labels/simplify things 
d0 = data.frame(cell_name= d0$cell_name, ribi.aucell = d0$ribi.aucell, rp.aucell=d0$rp.aucell, iesr.aucell=d0$iesr.aucell)
d30 = data.frame(cell_name= d30$cell_name, ribi.aucell = d30$ribi.aucell, rp.aucell=d30$rp.aucell, iesr.aucell=d30$iesr.aucell)



pbmc = readRDS("3004.naclt0.mostrecent.subtle.RDS")
pbmc@meta.data$cell_cycle = pbmc@active.ident
md.3004.t0 = pbmc@meta.data
md.3004.t0$cell_name = rownames(md.3004.t0)


#finalized t30 object where alphas have been removed
#pbmc = readRDS("3004.nacl.t30.20240527.RDS")
#md.3004.t30 = pbmc@meta.data
#md.3004.t30$cell_name = rownames(md.3004.t30)


#read in unfiltered data (from filtered bc matrix but not further filtered outside of cellranger)
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t0_2/NaCl_0.7M_3004_rep1_t0/filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
#dim(pbmc)
#6572 16160
unfiltered.3004.naclt0.md = pbmc@meta.data
unfiltered.3004.naclt0.md$cell_name = rownames(unfiltered.3004.naclt0.md)

#get indices of alphas
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))

#a_vec
indices.mfalpha = which(x = expr.mfalpha1 > 2)
#r_vec
indices.firstpass = which(x = expr.mfalpha1 <= 2)

#make a df of overall metadata before alphas and cluster 4 removed

df = pbmc@meta.data
df$expr.mfalpha1 = expr.mfalpha1$`MF%28ALPHA%291`
df$cell_name = rownames(df)
df$original.indices = 1:nrow(df)

#df with barcodes included in r_vec
df.sub = subset(df, df$expr.mfalpha1  <= 2)

dim(md.3004.t0)
dim(df.sub)
dd = subset(df.sub, !(df.sub$cell_name %in% md.3004.t0$cell_name))
dim(dd)

#dd is a df with barcodes corresponding to c_vec/represents the barcodes in the removed cluster

#use the unfiltered set of barcodes/indices and replace cell_cycle values with:
#dd$original.indices for the 'Unknown' 
#indices.mfalpha for the 'MATALPHA'




#first merge the unfiltered 3051 nacl t0 metadata with the aucell annotations from the seurat df.
#merge the assignments such that all existing ribi assignments go to the unfiltered pbmc metadata and N/As are added for 67 cells in seurat not in monocle 
result_merge.t0 <- merge(unfiltered.3004.naclt0.md, d0, by = "cell_name", all.x = TRUE)

#next, merge the seurat-based cell cycle labels from the df where alphas have already been removed 
#the cell_cycle column should contain discrete assignment from pca ananlysis in seurat as well as NAs for cells removed (here, just alphas)  
result_merge.test2 <- merge(result_merge.t0, md.3004.t0, by = "cell_name", all.x = TRUE)

#add named_data 
result_merge.test2$named_data = rep("NaCl_0.7M_3004_rep1_t0", length(rownames(result_merge.test2)))
result_merge.test2$n.umi = result_merge.test2$nCount_RNA.x
result_merge.test2$log2numi =  log(result_merge.test2$nCount_RNA.y, base = 2)

result_merge.test2$cell_cycle = as.character(result_merge.test2$cell_cycle)
result_merge.test2$cell_cycle[is.na(result_merge.test2$cell_cycle)] <- "MATALPHA"

#replace some of the 'matalpha' with the correct/'unknown' annotation 

#result_merge.test2$cell_cycle = replace()
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, dd$original.indices, "Unknown")

result_merge.test2$orig.ident.x <- NULL
result_merge.test2$nCount_RNA.x <- NULL
result_merge.test2$orig.ident.y <- NULL
result_merge.test2$nCount_RNA.y <- NULL
result_merge.test2$nFeature_RNA.y <- NULL
result_merge.test2$nCount_SCT <- NULL
result_merge.test2$nFeature_SCT <- NULL
result_merge.test2$nFeature_RNA.x  <- NULL


#cc labels are remade without spaaces or slashes
g2m.no.mating.indices = which(result_merge.test2$cell_cycle == "G2_M w/o mating")
g2m.mating.indices = which(result_merge.test2$cell_cycle == "G2_M + Mating")

mg1.mating.indices = which(result_merge.test2$cell_cycle == "M_G1 + Mating")
mg1.no.mating.indices = which(result_merge.test2$cell_cycle == "M_G1 w/o Mating")


result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, g2m.no.mating.indices, "G2_M_no_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, g2m.mating.indices, "G2_M_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, mg1.mating.indices, "M_G1_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, mg1.no.mating.indices, "M_G1_no_mating")

result_merge.test2.0 = result_merge.test2

write.table(result_merge.test2, file= "3004.nacl.statelabels.20240529.2.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)



########--------3004 nacl t30


#3004 nacl t0

#use this version, not the one without the (1)
df = read.delim("/Users/noahalexander/Downloads/NaCl_0.7M_3004_rep1.annotations (1).tsv")
d0 = subset(df, df$named_dataset =="NaCl_0.7M_3004_rep1_t0")
d30 = subset(df, df$named_dataset =="NaCl_0.7M_3004_rep1_t30")

#susbet to avoid retaining old cc labels/simplify things 
d0 = data.frame(cell_name= d0$cell_name, ribi.aucell = d0$ribi.aucell, rp.aucell=d0$rp.aucell, iesr.aucell=d0$iesr.aucell)
d30 = data.frame(cell_name= d30$cell_name, ribi.aucell = d30$ribi.aucell, rp.aucell=d30$rp.aucell, iesr.aucell=d30$iesr.aucell)



pbmc = readRDS("3004.naclt0.mostrecent.subtle.RDS")
pbmc@meta.data$cell_cycle = pbmc@active.ident
md.3004.t0 = pbmc@meta.data
md.3004.t0$cell_name = rownames(md.3004.t0)


#finalized t30 object where alphas have been removed
pbmc = readRDS("3004.nacl.t30.20240527.RDS")
pbmc@meta.data$cell_cycle = pbmc@active.ident
md.3004.t30 = pbmc@meta.data
md.3004.t30$cell_name = rownames(md.3004.t30)


#read in unfiltered data (from filtered bc matrix but not further filtered outside of cellranger)
pbmc_data <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t30_2/NaCl_0.7M_3004_rep1_t30/filtered_feature_bc_matrix/")
pbmc <- CreateSeuratObject(pbmc_data)
#dim(pbmc)
#6572 16160
unfiltered.3004.naclt30.md = pbmc@meta.data
unfiltered.3004.naclt30.md$cell_name = rownames(unfiltered.3004.naclt30.md)

#get indices of alphas
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))

#a_vec
indices.mfalpha = which(x = expr.mfalpha1 > 2)
#r_vec
indices.firstpass = which(x = expr.mfalpha1 <= 2)

#make a df of overall metadata before alphas and cluster 4 removed

df = pbmc@meta.data
df$expr.mfalpha1 = expr.mfalpha1$`MF%28ALPHA%291`
df$cell_name = rownames(df)
df$original.indices = 1:nrow(df)

#df with barcodes included in r_vec
df.sub = subset(df, df$expr.mfalpha1  <= 2)

dim(md.3004.t30)
dim(df.sub)
dd = subset(df.sub, !(df.sub$cell_name %in% md.3004.t30$cell_name))
dim(dd)

#dd is a df with barcodes corresponding to c_vec/represents the barcodes in the removed cluster

#use the unfiltered set of barcodes/indices and replace cell_cycle values with:
#dd$original.indices for the 'Unknown' 
#indices.mfalpha for the 'MATALPHA'




#first merge the unfiltered 3051 nacl t0 metadata with the aucell annotations from the seurat df.
#merge the assignments such that all existing ribi assignments go to the unfiltered pbmc metadata and N/As are added for 67 cells in seurat not in monocle 
result_merge.t30 <- merge(unfiltered.3004.naclt30.md, d30, by = "cell_name", all.x = TRUE)

#next, merge the seurat-based cell cycle labels from the df where alphas have already been removed 
#the cell_cycle column should contain discrete assignment from pca ananlysis in seurat as well as NAs for cells removed (here, just alphas)  
result_merge.test2 <- merge(result_merge.t30, md.3004.t30, by = "cell_name", all.x = TRUE)

#add named_data 
result_merge.test2$named_data = rep("NaCl_0.7M_3004_rep1_t30", length(rownames(result_merge.test2)))
result_merge.test2$n.umi = result_merge.test2$nCount_RNA.x
result_merge.test2$log2numi =  log(result_merge.test2$nCount_RNA.y, base = 2)

result_merge.test2$cell_cycle = as.character(result_merge.test2$cell_cycle)
result_merge.test2$cell_cycle[is.na(result_merge.test2$cell_cycle)] <- "MATALPHA"

#replace some of the 'matalpha' with the correct/'unknown' annotation 

#result_merge.test2$cell_cycle = replace()
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, dd$original.indices, "Unknown")

result_merge.test2$orig.ident.x <- NULL
result_merge.test2$nCount_RNA.x <- NULL
result_merge.test2$orig.ident.y <- NULL
result_merge.test2$nCount_RNA.y <- NULL
result_merge.test2$nFeature_RNA.y <- NULL
result_merge.test2$nCount_SCT <- NULL
result_merge.test2$nFeature_SCT <- NULL
result_merge.test2$nFeature_RNA.x  <- NULL


#cc labels are remade without spaaces or slashes
g2m.no.mating.indices = which(result_merge.test2$cell_cycle == "G2_M w/o mating")
g2m.mating.indices = which(result_merge.test2$cell_cycle == "G2_M + mating")

mg1.mating.indices = which(result_merge.test2$cell_cycle == "M_G1 + mating")
mg1.no.mating.indices = which(result_merge.test2$cell_cycle == "M_G1 w/o mating")


result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, g2m.no.mating.indices, "G2_M_no_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, g2m.mating.indices, "G2_M_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, mg1.mating.indices, "M_G1_mating")
result_merge.test2$cell_cycle = replace(result_merge.test2$cell_cycle, mg1.no.mating.indices, "M_G1_no_mating")

result_merge.test2.30 = result_merge.test2

#corresponding df for t0, just reading in a file for convenience 
t0.3004states = read.delim("3004.nacl.t0.statelabels.20240529.2.tsv")

test.out = rbind.fill( t0.3004states , result_merge.test2.30)

#test.out = rbind.fill( result_merge.test2.0 , result_merge.test2.30)


test.out$state = c(test.out$cell_cycle[1:18194], test.out$seurat_clusters[18195:nrow(test.out)])

saveRDS(test.out, file = "3004.NaCl.states.RDS")
write.table(test.out, file= "3004.NaCl.states.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)


write.table(result_merge.test2, file= "3004.nacl.t30.statelabels.20240529.2.tsv", quote = FALSE, col.names = T, sep = "\t", row.names = F)
