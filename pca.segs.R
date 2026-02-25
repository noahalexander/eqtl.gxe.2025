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


#######################plotting

p <- plot_cells(
  cds,
  reduction_method = "PCA",
  label_cell_groups = FALSE,
  color_cells_by = "sample"
)

p +
  labs(x = "PC1", y = "PC2") +  # Set axis labels
  theme(
    axis.title.x = element_text(size = 22),    # Increase x-axis label font size
    axis.title.y = element_text(size = 22),    # Increase y-axis label font size
    legend.title = element_text(size = 20),    # Increase legend title size
    legend.text = element_text(size = 18)      # Increase legend text/label size
  )


p <- plot_cells(
  cds,
  reduction_method = "PCA",
  label_cell_groups = FALSE,
  color_cells_by = "RiBi_Activity"
)

p +
  labs(x = "PC1", y = "PC2") +  # Set axis labels
  theme(
    axis.title.x = element_text(size = 22),    # Increase x-axis label font size
    axis.title.y = element_text(size = 22),    # Increase y-axis label font size
    legend.title = element_text(size = 20),    # Increase legend title size
    legend.text = element_text(size = 18)      # Increase legend text/label size
  )

p <- plot_cells(
  cds,
  reduction_method = "PCA",
  label_cell_groups = FALSE,
  color_cells_by = "iESR_Activity"
)

p +
  labs(x = "PC1", y = "PC2") +  # Set axis labels
  theme(
    axis.title.x = element_text(size = 22),    # Increase x-axis label font size
    axis.title.y = element_text(size = 22),    # Increase y-axis label font size
    legend.title = element_text(size = 20),    # Increase legend title size
    legend.text = element_text(size = 18)      # Increase legend text/label size
  )





