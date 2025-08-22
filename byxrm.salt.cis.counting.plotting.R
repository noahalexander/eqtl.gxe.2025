

df.nt0 = readRDS("repeat_fine_mapping/combined/A/nacl/t0/cis_only_test_CombinedResults.RDS")
df.nt0 = df.nt0$combined
#df = subset(df, df$FDR <= 0.05)
dim(df.nt0)

df.3 = readRDS("repeat_fine_mapping/combined/A/nacl/t30/cis_only_test_CombinedResults.RDS")
df.3 = df.3$combined
#df.3 = subset(df.3, df.3$FDR <= 0.05)
dim(df.3)


dd = merge(df, df.3, by = "transcript", all=T)
dim(dd)

df.spt0 = readRDS("repeat_fine_mapping/combined/A/sp/t0/cis_only_test_CombinedResults.RDS")
df.spt0= df.spt0$combined
#df.spt0 = subset(df.spt0, df.spt0$FDR <= 0.05)



df.t10 = readRDS("repeat_fine_mapping/combined/A/sp/t10/cis_only_test_CombinedResults.RDS")
df.t10 = df.t10$combined
#df.t10 = subset(df.t10, df.t10$FDR <= 0.05)


dl = merge(df.spt0, df.t10, by = "transcript", all=T)


dl = merge(df.nt0, df.3, by = "transcript", all=T)
dl = merge(df, df.3, by = "transcript", all=T)

dl = merge(df, df.t10, by = "transcript", all=T)
dl = merge(df, df.spt0, by = "transcript", all=T)


dl = merge(df.spt0, df.t10, by = "transcript", all=T)


###########-----------------------------

dl = merge(df.nt0, df.3, by = "transcript", all=T)

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
    vec[i] = "t30.only"
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
  scale_color_manual(values = c("both" = "red", "t30.only" = "blue", "t0.only" = "green", "neither" = "black"))



#add de stats from seurat to transcripts with a sig pval and add na to rest 


timepoint.de.markers$transcript = rownames(timepoint.de.markers)
dl2 = merge(dl, timepoint.de.markers, by ="transcript", all=T)

ggplot(data = dl2, aes(log(LOD.x), log(LOD.y))) +
  geom_point(aes(color = avg_log2FC)) +
  scale_color_gradient(low = "darkblue", high = "red") 

dd = subset(dl2, dl2)








################3


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

#####


dl = merge(df, df.spt0, by = "transcript", all=T)

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
    vec[i] = "sp.t0.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <=0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "nacl.t0.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] > 0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "neither"
  }
}

dl$vec = vec
ggplot(data = dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point() +
  scale_color_manual(values = c("both" = "red", "nacl.t0.only" = "blue", "sp.t0.only" = "green", "neither" = "black"))



#merging the de df from seurat
#################3051 nacl t0

pbmc_data.sp.t0 <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/seurat_in")
pbmc.sp.1 <- CreateSeuratObject(pbmc_data.sp.t0)



#################3051 sp t10

pbmc_data.sp.t10 <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t10//filtered_feature_bc_matrix/seurat_in")
pbmc.sp.2 <- CreateSeuratObject(pbmc_data.sp.t10)







#################3051 nacl t0

pbmc_data.1 <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t0///filtered_feature_bc_matrix/")
pbmc.1 <- CreateSeuratObject(pbmc_data.1)
dim(pbmc.1)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc.1@meta.data$nCount_RNA, base = 2))



#################3051 nacl t30
pbmc_data.2 <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t30/outs//filtered_feature_bc_matrix/")
pbmc.2 <- CreateSeuratObject(pbmc_data.2)
dim(pbmc.2)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc.2@meta.data$nCount_RNA, base = 2))




pbmc <- merge(pbmc.1, y = pbmc.2, add.cell.ids = c("t0", "t30"), project = "nacl_byxrm")


pbmc = pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)







dd = read.table("/Users/noahalexander/Downloads/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/features.tsv.gz")



RNA <- pbmc@assays$RNA
RNA@counts@Dimnames[[1]] = dd$V1
RNA@data@Dimnames[[1]] = dd$V1

rownames(RNA@meta.features) =dd$V1
pbmc@assays$RNA <- RNA





tvec = c(rep("t0",16160), rep("t30", 14664))
pbmc@meta.data$timepoint = tvec

Idents(pbmc) <- "timepoint"

timepoint.de.markers <- FindMarkers(pbmc, ident.1 = "t0", ident.2 = "t30")
