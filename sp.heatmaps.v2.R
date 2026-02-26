library(Seurat)
library(sctransform)
library(dplyr)
library(stringr)
library(ggplot2)

df = read.csv("cc.genes.df.csv")
df2 = read.delim("pbio.2004050.csv", sep = ",")
esr.genes = word(df2$Annotation, 4)



pbmc_data <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t0/filtered_feature_bc_matrix/seurat_in")
pbmc <- CreateSeuratObject(pbmc_data)
dim(pbmc)
#pbmc <- NormalizeData(pbmc)
#dim(pbmc)
hist(log(pbmc@meta.data$nCount_RNA, base = 2))

#do things change if i do this after sctransform or does it use same slot regardless?
expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)

#expr.mfalpha2<- FetchData(object = pbmc, vars = c("MF(ALPHA)2"))
dim(pbmc)
pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
dim(pbmc)
#pbmc = pbmc[, which(x = expr.mfalpha2 < .01)]
#dim(pbmc)


pbmc= SCTransform(pbmc)
dim(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
dim(pbmc)
#unlike above, do not use just cc genes for pca
pbmc= RunPCA(pbmc, npcs = 12)
dim(pbmc)
pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca")
dim(pbmc)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.25  )
dim(pbmc)


FeaturePlot(pbmc, features = cc_markers)
DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)


#current
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

dim(top10)


#new
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  slice_head(n = 15) %>%
  ungroup() -> top10

dim(top10)

DoHeatmap(pbmc, features = c(top10$gene, cc_markers, "SPG1")) 


