library(Seurat)
library(sctransform)
library(dplyr)
library(stringr)
library(ggplot2)

df = read.csv("cc.genes.df.csv")
df2 = read.delim("pbio.2004050.csv", sep = ",")
esr.genes = word(df2$Annotation, 4)

######sp t0

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


#####################complexheatmap version 

ht_opt$message = FALSE
DefaultAssay(pbmc) <- "SCT"

features <- unique(c(top10$gene, cc_markers, "SPG1"))

scaled_genes <- rownames(pbmc[["SCT"]]@scale.data)
if (length(scaled_genes) == 0 || any(!features %in% scaled_genes)) {
    pbmc <- ScaleData(pbmc, assay = "SCT", features = features, verbose = FALSE)
}

mat_all <- GetAssayData(pbmc, assay = "SCT", slot = "scale.data")
features_use <- intersect(features, rownames(mat_all))
mat <- as.matrix(mat_all[features_use, , drop = FALSE])

clusters <- Idents(pbmc)
ord <- order(as.integer(clusters))
mat <- mat[, ord, drop = FALSE]
clusters <- droplevels(clusters[ord])

roman_labels <- c("I","II","III","IV","V","VI")
stopifnot(length(levels(clusters)) <= length(roman_labels))
levels(clusters) <- roman_labels[seq_along(levels(clusters))]

lims <- c(-1.5, 1.5)
mat_clip <- pmax(pmin(mat, lims[2]), lims[1])

p <- 0.6
x <- mat_clip / lims[2]
mat_plot <- sign(x) * (abs(x)^p) * lims[2]

col_fun <- colorRamp2(
    c(lims[1], 0, lims[2]),
    c("white", "#FFE066", "#1A0F0B")
)

ht <- Heatmap(
    mat_plot,
    name = "Expression",
    col = col_fun,
    column_split = clusters,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 13),
    column_title_gp = gpar(fontsize = 16, fontface = "bold"),
    use_raster = FALSE,
    heatmap_legend_param = list(
        at = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
        
        # Use stable title position
        title_position = "topleft",
        
        # Slightly smaller so it fits RStudio device
        title_gp = gpar(fontsize = 16, fontface = "bold"),
        labels_gp = gpar(fontsize = 14),
        
        legend_height = unit(80, "mm"),
        legend_width  = unit(20, "mm")
    )
)

grid.newpage()
draw(ht, heatmap_legend_side = "right")


