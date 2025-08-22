df.spt0 = readRDS("repeat_fine_mapping/combined/A/sp/t0/cis_only_test_CombinedResults.RDS")
df.spt0= df.spt0$combined
#df.spt0 = subset(df.spt0, df.spt0$FDR <= 0.05)



df.t10 = readRDS("repeat_fine_mapping/combined/A/sp/t10/cis_only_test_CombinedResults.RDS")
df.t10 = df.t10$combined


dl = merge(df.spt0, df.t10, by = "transcript", all=T)

dl$FDR.x[is.na(dl$FDR.x)] <- 1
dl$FDR.y[is.na(dl$FDR.y)] <- 1

vec= vector()
for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <= 0.05 & dl$FDR.y[i] <= 0.05) {
    vec[i] = "both"
  }
}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] >0.05 & dl$FDR.y[i] <= 0.05) {
    vec[i] = "t10.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <=0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "t0.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] > 0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "neither"
  }
}

dl$vec = vec
ggplot(data = dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point() +
  scale_color_manual(values = c("both" = "red", "t10.only" = "blue", "t0.only" = "green", "neither" = "black"))

#############

pbmc_data.sp.t0 <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/seurat_in")
pbmc.sp.1 <- CreateSeuratObject(pbmc_data.sp.t0)



#################3051 sp t10

pbmc_data.sp.t10 <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t10//filtered_feature_bc_matrix/seurat_in")
pbmc.sp.2 <- CreateSeuratObject(pbmc_data.sp.t10)


pbmc <- merge(pbmc.sp.1, y = pbmc.sp.2, add.cell.ids = c("t0", "t10"), project = "nacl_sp")


pbmc = pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)





dd = read.table("/Users/noahalexander/Downloads/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/features.tsv.gz")



RNA <- pbmc@assays$RNA
RNA@counts@Dimnames[[1]] = dd$V1
RNA@data@Dimnames[[1]] = dd$V1

rownames(RNA@meta.features) =dd$V1
pbmc@assays$RNA <- RNA



tvec = c(rep("t0",10472), rep("t10",19682 ))
pbmc@meta.data$timepoint = tvec

Idents(pbmc) <- "timepoint"

timepoint.de.markers.sp.t10vst0 <- FindMarkers(pbmc, ident.1 = "t0", ident.2 = "t10")




timepoint.de.markers.sp.t10vst0$transcript = rownames(timepoint.de.markers.sp.t10vst0)
dl2 = merge(dl, timepoint.de.markers.sp.t10vst0, by ="transcript", all=T)

ggplot(data = dl2, aes(log(LOD.x), log(LOD.y))) +
  geom_point(aes(color = avg_log2FC)) +
  scale_color_gradient(low = "darkblue", high = "red") 
