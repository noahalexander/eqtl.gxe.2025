df.0 = readRDS("repeat_fine_mapping/combined/3004/nacl/t0/cis_only_test_CombinedResults.RDS")
df.0= df.0$combined
#df.spt0 = subset(df.spt0, df.spt0$FDR <= 0.05)



df.30 = readRDS("repeat_fine_mapping/combined/3004//nacl/t30/cis_only_test_CombinedResults.RDS")
df.30 = df.30$combined


dl = merge(df.0, df.30, by = "transcript", all=T)

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
    vec[i] = "nacl.t0.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <=0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "nacl.t30.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] > 0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "neither"
  }
}

dl$vec = vec
ggplot(data = dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point() +
  scale_color_manual(values = c("both" = "red", "nacl.t0.only" = "blue", "nacl.t30.only" = "green", "neither" = "black"))


#################cbsxyjm nacl t0

pbmc_data.1 <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t0_2/NaCl_0.7M_3004_rep1_t0/filtered_feature_bc_matrix/")
pbmc.1 <- CreateSeuratObject(pbmc_data.1)
dim(pbmc.1)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc.1@meta.data$nCount_RNA, base = 2))



#################cbsxyjm nacl t30
pbmc_data.2 <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t30_2/NaCl_0.7M_3004_rep1_t30/filtered_feature_bc_matrix/")
pbmc.2 <- CreateSeuratObject(pbmc_data.2)
dim(pbmc.2)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc.2@meta.data$nCount_RNA, base = 2))


pbmc <- merge(pbmc.1, y = pbmc.2, add.cell.ids = c("t0", "t30"), project = "nacl")


pbmc = pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)







dd = read.table("/Users/noahalexander/Downloads/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/features.tsv.gz")



RNA <- pbmc@assays$RNA
RNA@counts@Dimnames[[1]] = dd$V1
RNA@data@Dimnames[[1]] = dd$V1

rownames(RNA@meta.features) =dd$V1
pbmc@assays$RNA <- RNA





tvec = c(rep("t0",18194), rep("t30",18026 ))
pbmc@meta.data$timepoint = tvec

Idents(pbmc) <- "timepoint"

timepoint.de.markers.3004 <- FindMarkers(pbmc, ident.1 = "t0", ident.2 = "t30")



timepoint.de.markers.3004$transcript = rownames(timepoint.de.markers.3004)
dl2 = merge(dl, timepoint.de.markers.3004, by ="transcript", all=T)

ggplot(data = dl2, aes(log(LOD.x), log(LOD.y))) +
  geom_point(aes(color = avg_log2FC)) +
  scale_color_gradient(low = "darkblue", high = "red") 
