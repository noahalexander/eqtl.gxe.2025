library(monocle3)
library(ggplot2)
library(viridis)
library(patchwork)

cds.3051.nacl.t0 <- load_mm_data(mat_path = "/Users/noahalexander/seg_cellranger_report/t0/filtered_feature_bc_matrix/matrix.mtx", 
  feature_anno_path = "/Users/noahalexander/seg_cellranger_report/t0/filtered_feature_bc_matrix/features.tsv", 
  cell_anno_path = "/Users/noahalexander/seg_cellranger_report/t0/filtered_feature_bc_matrix/barcodes.tsv")


cds.3051.nacl.t30 <- load_mm_data(mat_path = "/Users/noahalexander/seg_cellranger_report/t30/filtered_feature_bc_matrix/matrix.mtx", 
  feature_anno_path = "/Users/noahalexander/seg_cellranger_report/t30/filtered_feature_bc_matrix/features.tsv", 
  cell_anno_path = "/Users/noahalexander/seg_cellranger_report/t30/filtered_feature_bc_matrix/barcodes.tsv")



cds.3051.sp.t0 <- load_mm_data(mat_path = "/Users/noahalexander/q.seg.2023/t0/filtered_feature_bc_matrix/matrix.mtx", 
  feature_anno_path = "/Users/noahalexander/q.seg.2023/t0/filtered_feature_bc_matrix/features.tsv", 
  cell_anno_path = "/Users/noahalexander/q.seg.2023/t0/filtered_feature_bc_matrix/barcodes.tsv")


cds.3051.sp.t10 <- load_mm_data(mat_path = "/Users/noahalexander/q.seg.2023/t10/filtered_feature_bc_matrix/matrix.mtx", 
  feature_anno_path = "/Users/noahalexander/q.seg.2023/t10/filtered_feature_bc_matrix/features.tsv", 
  cell_anno_path = "/Users/noahalexander/q.seg.2023/t10/filtered_feature_bc_matrix/barcodes.tsv")



cds = combine_cds(list(cds.3051.nacl.t0, cds.3051.nacl.t30, cds.3051.sp.t0, cds.3051.sp.t10))
cds@colData$sample = as.factor(cds@colData$sample)
cds@colData$sample.og = as.factor(cds@colData$sample)


colData(cds)$sample <- dplyr::recode(
  as.character(colData(cds)$sample),
  '1' = 'NaCl t0',
  '2' = 'NaCl t30',
  '3' = 'Stationary Phase t0',
  '4' = 'Stationary Phase t10'
  # Add more mappings as needed
)


cds = preprocess_cds(cds, use_genes = intersect(rownames(cds), botstein.naming$ORF))

cds = mono.aucell(cds)

plot_cells(cds, reduction_method = "PCA", label_cell_groups = F, color_cells_by = "sample")

cds@colData$Sample = cds@colData$sample
cds@colData$RiBi_Activity = cds@colData$ribi.aucell
cds@colData$RP_Activity = cds@colData$rp.aucell
cds@colData$iESR_Activity = cds@colData$iesr.aucell


#######################plotting

#samples
p <- plot_cells(
  cds,
  reduction_method = "PCA",
  label_cell_groups = FALSE,
  color_cells_by = "Sample"
)

p +
  labs(x = "PC1", y = "PC2") +  # Set axis labels
  theme(
    axis.title.x = element_text(size = 22),    # Increase x-axis label font size
    axis.title.y = element_text(size = 22),    # Increase y-axis label font size
    legend.title = element_text(size = 20),    # Increase legend title size
    legend.text = element_text(size = 18)      # Increase legend text/label size
  )

#ribi aucell values
p <- plot_cells(
    cds,
    reduction_method = "PCA",
    label_cell_groups = FALSE,
    color_cells_by = "RiBi_Activity"
)

p +
    labs(x = "PC1", y = "PC2") +
    guides(colour = guide_colourbar(title = "RiBi Activity")) +  # <- change title only
    theme(
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)
    )

#iesr aucell values

p <- plot_cells(
  cds,
  reduction_method = "PCA",
  label_cell_groups = FALSE,
  color_cells_by = "iESR_Activity"
)

p +
  labs(x = "PC1", y = "PC2") +
  guides(colour = guide_colourbar(title = "iESR Activity")) +  # <- change title only
  theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )


#rp aucell values

p <- plot_cells(
    cds,
    reduction_method = "PCA",
    label_cell_groups = FALSE,
    color_cells_by = "RP_Activity"
)

p +
    labs(x = "PC1", y = "PC2") +
    guides(colour = guide_colourbar(title = "RP Activity")) +  # <- change title only
    theme(
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)
    )

#######################

theme_big <- theme(
  axis.title.x = element_text(size = 22),
  axis.title.y = element_text(size = 22),
  legend.title = element_text(size = 20),
  legend.text  = element_text(size = 18),
  plot.margin  = margin(5.5, 5.5, 5.5, 5.5)
)

# Make your four plots (as you already do)
p_sample <- plot_cells(cds, reduction_method="PCA", label_cell_groups=FALSE, color_cells_by="Sample") +
  labs(x="PC1", y="PC2") + theme_big

p_ribi <- plot_cells(cds, reduction_method="PCA", label_cell_groups=FALSE, color_cells_by="RiBi_Activity") +
  labs(x="PC1", y="PC2") +
  guides(colour = guide_colourbar(title = "RiBi Activity")) +
  theme_big

p_iesr <- plot_cells(cds, reduction_method="PCA", label_cell_groups=FALSE, color_cells_by="iESR_Activity") +
  labs(x="PC1", y="PC2") +
  guides(colour = guide_colourbar(title = "iESR Activity")) +
  theme_big

p_rp <- plot_cells(cds, reduction_method="PCA", label_cell_groups=FALSE, color_cells_by="RP_Activity") +
  labs(x="PC1", y="PC2") +
  guides(colour = guide_colourbar(title = "RP Activity")) +
  theme_big

# Helper: lock a fixed legend column width so all panels match
panel_with_legend <- function(p, legend_rel_width = 0.30) {
  lg <- get_legend(p + theme(legend.position = "right"))
  p_nolg <- p + theme(legend.position = "none")
  plot_grid(p_nolg, lg, nrow = 1, rel_widths = c(1, legend_rel_width))
}

# Wrap each panel with a legend of the SAME relative width
w_sample <- panel_with_legend(p_sample, legend_rel_width = 0.35)  # discrete legend often needs a bit more
w_ribi   <- panel_with_legend(p_ribi,   legend_rel_width = 0.35)
w_iesr   <- panel_with_legend(p_iesr,   legend_rel_width = 0.35)
w_rp     <- panel_with_legend(p_rp,     legend_rel_width = 0.35)

# Now assemble the 2x2 grid with aligned panels
fig4 <- plot_grid(
    w_sample, w_ribi,
    w_iesr,   w_rp,
    ncol = 2,
    align = "hv",
    axis = "tblr",
    labels = c("A", "B", "C", "D"),   # or c("a","b","c","d")
    label_size = 18,
    label_fontface = "bold"
)

fig4


ggsave(
  "Figure_PCA_4panel.png",
  fig4,
  width = 14, height = 10, dpi = 300,
  bg = "white"
)


