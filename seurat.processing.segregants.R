library(Seurat)
library(sctransform)
library(dplyr)
library(stringr)

df = read.csv("cc.genes.df.csv")
df2 = read.delim("pbio.2004050.csv", sep = ",")
esr.genes = word(df2$Annotation, 4)

#new.markers = c("CLB4","WHI3","CIS3","CDC28")
#more.markers = c("SIC1","CLB2", "CLN3", "SWI4", "CDC6", "CDC47", "KAR4", "AGA2", "SST2", "FUS1", "FAR1", "DBF2", "KIN3", "TEC1", "PCL9", "SST2", "CDC2", "POL1", "POL2", "MSH2", "MCM2", "MCM3", "DUB4", "GIC1", "CMN67", "NUF1","CIN8", "KAR9","TUB1", "TUB3")
#clb2.cluster = 
#more.markers = c("SIC1","CLB2", "CLN3", "SWI4", "CDC6", "CDC47", "KAR4", "AGA2", "SST2", "FUS1", "FAR1", "DBF2", "KIN3", "TEC1", "PCL9", "SST2", "CDC2", "POL1", "POL2", "MSH2", "MCM2", "MCM3", "DUB4", "GIC1", "CMN67", "NUF1","CIN8", "KAR9","TUB1", "TUB3")
#more.markers = c("CLB2", "CLN3", "SWI4", "CDC6", "CDC47", "KAR4", "AGA2", "SST2", "FUS1", "FAR1", "DBF2", "KIN3", "TEC1", "PCL9", "SST2", "CDC2", "POL1", "POL2", "MSH2", "MCM2", "MCM3", "DUB4", "GIC1", "CMN67", "NUF1","CIN8", "KAR9","TUB1", "TUB3", "ARG1", "ARG3", "ARG5", "ARG8")

#all.marker = c(more.markers, "MET2", "ECM1", "CLB2", 'CDC5', "CDC20",  "SWI5", "WSC4", "PMP1", "PMA1" , "PMA2", "CLB2", "CLB1", "CLB2", "SWI5",  'BUD4',"ALG7", "FKS1", 'GAS1', 'GOG5', 'PMT1',  'PMI40', "ARG1", "ARG3", "ARG5", "ARG8", "MET17", "MET2", "ECM17", "CLB2", "CDC5", "CDC20",  "SWI5", "WSC4", "PMP1","PMA1" , "PMA2",  "CLB2",
#  "CLB1", "CLB2", "SWI5", 'BUD4', 'ASE1', "KIP2", "MOB1", "NUM1", 'YCL012W', "BUD3", "CHA1", "YCL063W", "YLR057W",  "YML033W" , "MCM2", 'MCM3', "CDC54", "CDC46", 'MCM6',  "CDC47", 'KAR4', "AGA1", "SST2",  "FUS1", "CLB4", "BUD3", "CPR8", "PRO2", "YCL012W", "YCL063W", "YGL217C", "YNL043C", "YDR130C",  'YOL030W', "CDC2", "POL1",  "POL2", "CDC21", "CDC45", "PMS1",
#  "MSH2", "MCM2", "MCM3", "CDC54", "CDC47", "MCM6", "CDC47",  "CDC54",  "CDC6", "BUD3", "BUD4", "BUD8", 'BUD9', "BEM1", 'GIC1', 'MSB1', 'MSB2',
#"CNM67", "NUF1", "SPC42", "SPC97", 'SPC98', 'TUB4', "SPC34", "BIM1", "BUB1", "IPL1", 'KAR3', "SLK19", "CIN8", "KAR9", "KIP1", "STU2",  "VIK1", 'BUB2', "CIK1", 'KIP2', "KIP3",  "NUM1", "TUB2", 'ASE1', "AUA1",  "MEP3", "HXT1", "RGT2", "FET3" , "FTR1", "PHO3" , "PHO8", "PS4",  "SSP2")


#MET2. Finally, ECM1
#CLB2, CDC5, CDC20, and SWI5. 
#WSC4, PMP1
#PMA1 and PMA2. The CLB2
#CLB1, CLB2, SWI5, and BUD4,
#ALG7, FKS1, GAS1, GOG5, PMT1, and PMI40,
#
#ARG1, ARG3, ARG5,6 and ARG8)
#more.markers = c("CLB2", "CLN3", "SWI4", "CDC6", "CDC47", "KAR4", "AGA2", "SST2", "FUS1", "FAR1", "DBF2", "KIN3", "TEC1", "PCL9", "SST2", "CDC2", "POL1", "POL2", "MSH2", "MCM2", "MCM3", "DUB4", "GIC1", "CMN67", "NUF1","CIN8", "KAR9","TUB1", "TUB3", "ARG1", "ARG3", "ARG5", "ARG8")
#ALG7, FKS1, GAS1, GOG5, PMT1, and PMI40,
#MET17, MET2, ECM17
#CLB2, CDC5, CDC20, and SWI5.
#WSC4, PMP1,PMA1 and PMA2. The CLB2
#CLB1, CLB2, SWI5, and BUD4
#ASE1
#KIP2, MOB1, NUM1, YCL012W, BUD3, CHA1, YCL063W, YLR057W, and YML033W
#MCM2, MCM3, CDC54, CDC46, MCM6, and CDC47
#KAR4, AGA1, SST2, and FUS1
#CLB4, BUD3, CPR8, PRO2, YCL012W, YCL063W, YGL217C, YNL043C, YDR130C, and YOL030W;
#CDC2, POL1, and POL2
#CDC21), and genes involved in initiation of DNA synthesis (e.g. CDC45) PMS1 and MSH2 
#(MCM2, MCM3, CDC54, CDC47, MCM6, CDC47, and CDC54) and CDC6
#BUD3, BUD4, BUD8, BUD9, BEM1, GIC1, MSB1, and MSB2)
#CNM67, NUF1, SPC42, SPC97, SPC98, and TUB4), one (SPC34) 
#BIM1, BUB1, IPL1, KAR3, and SLK19#
#(CIN8, KAR9, KIP1, STU2, and VIK1) peak. Five genes peak during G2 (BUB2, CIK1, KIP2, KIP3, and NUM1), as well as the major 􏰂 tubulin TUB2. Finally, one gene (ASE1) reaches peak expres- sion during M.
#(AUA1 and MEP3), sugars (e.g., HXT1 and RGT2), and iron (FET3 and FTR1). We also identified the acid phosphatases (e.g., PHO3 and PHO8).
#PS4 and SSP2,

#Nearly all of these genes reach peak expression late in the cell cycle during M and M/G1.
#(GAP1), ammonia (AUA1 and MEP3), sugars (e.g., HXT1 and RGT2), and iron (FET3 and FTR1). We also identified the acid phosphatases (e.g., PHO3 and PHO8
  
  #A number of genes associated with functions in spe- cialized developmental pathways show cell cycle reg- ulation. These include the apparently sporulation-spe- cific genes
#  SPS4 and SSP2, which have peak expression in M and M/G1, respectively. These might represent cases such as SPO12,
  
#  SMC3, MCD1, PDS1, and PDS5
#  PMA1
  
#  mat cluster?
#    Y
#  fks1
#  histone
#  met
#  clb2
#  mcm
#  sic1
  
  
  
  
  
  
  
  ####3051 t0
  pbmc_data <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t0/filtered_feature_bc_matrix/")
  pbmc <- CreateSeuratObject(pbmc_data)
  dim(pbmc)
  #pbmc <- NormalizeData(pbmc)
  #dim(pbmc)
  
  
  expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
  expr.dse2 <- FetchData(object = pbmc, vars = c("DSE2"))
  #expr.mfalpha2<- FetchData(object = pbmc, vars = c("MFA2"))
  hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)
  
  dim(pbmc)
  pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
  dim(pbmc)
  #pbmc = pbmc[, which(x = expr.mfalpha2 < 1)]
  dim(pbmc)
  
  hist(log(pbmc@meta.data$nCount_RNA, base = 2))
  
  pbmc= SCTransform(pbmc) 
  dim(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
  dim(pbmc)
  #pca using cc genes
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12) 
  dim(pbmc)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3 )
  dim(pbmc)
  
  ElbowPlot(pbmc)
  
  FeaturePlot(pbmc, features = cc_markers)
  FeaturePlot(pbmc, features = cc_markers, dims = c(2,3))
  
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10, dims = c(2,3))
  
  #mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
  
  #######
  #new.cluster.ids <- c()
  #names(new.cluster.ids) <- levels(pbmc)
  #pbmc <- RenameIdents(pbmc, new.cluster.ids)
  
  
  #######################
  
  
  
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() -> top10
  
  
  #making pheatmap heatmaps instead of seurat ones where cols can't be clustered
  #df = pbmc@assays$SCT@scale.data
  #df = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")  
  #df = subset(df, rownames(df) %in% top10$gene)
  #dd = pbmc@meta.data
  #dd = data.frame(seurat_clusters = dd$seurat_clusters, row.names = rownames(dd))
  #pheatmap(df, show_colnames = F, annotation_col = dd)
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers)) 
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers, more.markers)) 
  
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers, all.marker)) 
  
  
  #pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  #pbmc.markers %>%
  #  group_by(cluster) %>%
  #  dplyr::filter(avg_log2FC > 1.5) %>%
  #  slice_head(n = 5) %>%
  #  ungroup() -> top10
  #FeaturePlot(pbmc, features = top10$gene)
  
  # pbmc = readRDS("3051.seurat.nacl.t0.20240417.RDS")
  
  
  
  
  #naclt03051.cc = c("G2/M", "S", "G1/S", "S", "M/G1","G2/M",  "M/G1", "G1" )
  #naclt03051.cc.simple = c("G2/M + Mating", "S/G2", "M/G1 + Mating", "G2/M no Mating", "S", "G1/S", "M/G1 mating", "M/G1 no Mating" )
  naclt03051.cc = c("G2_M + mating", "S", "G1_S", "S", "G1","G2_M w/o mating",  "M_G1 + mating", "M/G1 w/o mating")
  

  names(naclt03051.cc) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, naclt03051.cc)
  
#DoHeatmap(pbmc, features = c(top10$gene, cc_markers), group.by = c("M_G1 + mating", "M/G1 w/o mating", "G1", "G1_S", "S", "G2_M w/o mating","G2_M + mating" )) 
  
  
  #make new pca plots using collapsed labels 
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10, dims = c(2,3))
  
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers)) + scale_fill_gradientn(colors = c( "white", "gray", "black"))
  DoHeatmap(object = pbmc, features = c(top10$gene, cc_markers)) + scale_fill_gradientn(colors = c( "white", "lemonchiffon", "wheat4", "black"))
  
  
  levels(pbmc)
  #intially: "G2_M + mating"   "S"               "G1_S"            "G1"              "G2_M w/o mating" "M_G1 + mating"   "M/G1 w/o mating"
  
  levels(pbmc) = c("M_G1 + mating", "M/G1 w/o mating", "G1", "G1_S", "S", "G2_M w/o mating","G2_M + mating" )
  levels(pbmc)
  
  
  #generate 'final' heatmap with ordered states/clusters 
  DoHeatmap(object = pbmc, features = c(top10$gene, cc_markers)) + scale_fill_gradientn(colors = c( "white", "lemonchiffon", "wheat4", "black"))
  

 
  
  
  ##########################3051 nacl t30
  pbmc_data <- Read10X(data.dir = "/Users/noahalexander/3051.correctref.nacl.segregants/t30/outs//filtered_feature_bc_matrix/")
  pbmc <- CreateSeuratObject(pbmc_data)
  dim(pbmc)
  #pbmc <- NormalizeData(pbmc)
  #dim(pbmc)
  hist(log(pbmc@meta.data$nCount_RNA, base = 2))
  
  
  #do things change if i do this after sctransform or does it use same slot regardless?
  expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
  hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)
  
  #expr.mfalpha2<- FetchData(object = pbmc, vars = c(""))
  dim(pbmc)
  pbmc = pbmc[, which(x = expr.mfalpha1 <= 2)]
  dim(pbmc)
  #pbmc = pbmc[, which(x = expr.mfalpha2 < .01)]
  #dim(pbmc)
  
  
  pbmc= SCTransform(pbmc) 
  dim(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
  dim(pbmc)
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12) 
  dim(pbmc)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.55 )
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.43 )
  
  
  
  #pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5 )
  
  dim(pbmc)
  
  
  FeaturePlot(pbmc, features = cc_markers)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  
  # pbmc = readRDS("3051.nacl.t30.seurat.RDS")
  
  #mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
  #naclt303051.cc = c("G2/M", "G1/S", "G1", "M/G1", "S","G2/M",  "G2/M" )
  
  #naclt303051.cc = c("G2/M", "G1/S", "G1", "M/G1", "S","G2/M",  "G2/M" )
  
  #for 0.5 res clusters in 12 dim pc space - crude version
  #naclt303051.cc = c("M/G1", "G1/S", "G2/M", "M/G1", "G2/M","S",  "G2/M" , "G1/S", "Stress","S")
  
  
  #for 0.5 res clusters in 12 dim pc space - more nuanced version
  #naclt303051.cc = c("light M/G1 with stress", "G1/S", "G1", "stronger M/G1", "G2/M","S",  "G2/M" , "G1/S", "Stress","S with strong RP")
  
  #naclt303051.cc = c("G1", "G1/S", "G2/M", "M/G1", "G2/M","S",  "G2/M" , "G1/S", "Stress","S + RP")
  #naclt303051.cc = c("G1", "G1/S", "G2/M", "M/G1", "G2/M","S",  "G2/M" , "G1/S", "Stress","S")
  
  #for 0.55 res clustering in cc pc space with 12 dim
  
  
  #read in "3051.naclt30.res0.55.RDS"
  #for res 0.55 
  #naclt303051.cc = c("G2_M", "G1", "G1_S", "M_G1", "S","G2_M",  "G2_M" , "G1" ,"Stress","S")
  #for res 0.43
  #naclt303051.cc = c("G1_S", "G1", "G2_M",  "G2_M" , "M_G1" ,"S", "G2_M","S")
  
  #for 0.43 res
  naclt303051.cc = c("G1_S", "G1", "G2_M_no_mating",  "G2_M_mating" , "M_G1" ,"S", "G2_M_mating","S")
  
  
  
  
  names(naclt303051.cc) <- levels(pbmc)
  
  pbmc <- RenameIdents(pbmc, naclt303051.cc)
  
  
  
  #######################
  
  
  #new markers for collapsed clusters
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top10
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers, "MFA2", "CTT1"))
  DoHeatmap(pbmc, features = cc_markers)
  
  #3051.nacl.t30.seurat.RDS
  
  
  #################3004 nacl t0 
  pbmc_data <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t0_2/NaCl_0.7M_3004_rep1_t0/filtered_feature_bc_matrix/")
  pbmc <- CreateSeuratObject(pbmc_data)
  dim(pbmc)
  #pbmc <- NormalizeData(pbmc)
  #dim(pbmc)
  
  #do things change if i do this after sctransform or does it use same slot regardless?
  expr.mfalpha1 <- FetchData(object = pbmc, vars = c("MF%28ALPHA%291"))
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
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
  dim(pbmc)
  #pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3,  )
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3,  )
  
  
  dim(pbmc)
  
  hist(log(pbmc@meta.data$nCount_RNA, base = 2))
  hist(log(expr.mfalpha1$`MF%28ALPHA%291`, base = 2), breaks = 20)
  
  FeaturePlot(pbmc, features = cc_markers)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  
  #mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top10
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers))
  
  
  
  
  
  #pbmc = readRDS("3051.nacl.t30.seurat.RDS")
  
  #removal of cluster of indistinct cells that overlap with 'unclear' ones 
  dim(pbmc)
  cells_to_remove <- WhichCells(pbmc, idents = c(4))
  pbmc <- subset(pbmc, cells = setdiff(Cells(pbmc), cells_to_remove))
  dim(pbmc)
  
  
  #repeat sctransform using subset object
  pbmc= SCTransform(pbmc) 
  dim(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
  dim(pbmc)
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
  dim(pbmc)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3,  )
  dim(pbmc)
  
  
  
  FeaturePlot(pbmc, features = cc_markers)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  
  
  #######################
  #find markers and make another heatmap using final (cc) clusters
  #pbmc = readRDS("3004.nacl.t0.alphaandcluster4filtered.RDS")
  
  
  
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top10
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers))
  
  
  #new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
  #names(new.cluster.ids) <- levels(pbmc)
  #pbmc <- RenameIdents(pbmc, new.cluster.ids)
  #DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  
  #naclt03004.cc = c("G2/M", "S", "M/G1", "G2/M", "S", "G1/S", "G1" )
  #naclt03004.cc.simple = c("G2/M", "S/G2", "M/G1 + Mating", "G2/M no Mating", "S", "G1/S", "S/G2", "M/G1 no Mating" )
  naclt03004.cc = c("G2_M + mating", "S", "M_G1 + Mating", "G2_M w/o mating", "S", "G1_S", "S", "M_G1 w/o Mating" )
  
  
  
  
  names(naclt03004.cc) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, naclt03004.cc)
  
  
  #########################
  #########################
  pbmc_data <- Read10X(data.dir = "/Users/noahalexander/NaCl_0.7M_3004_rep1/NaCl_0.7M_3004_rep1_t30_2/NaCl_0.7M_3004_rep1_t30//filtered_feature_bc_matrix/")
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
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
  dim(pbmc)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5,  )
  dim(pbmc)
  
  
  FeaturePlot(pbmc, features = cc_markers)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10, dims = c(2,3))
  
  
  #mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top10
  
  DoHeatmap(pbmc, features = c(c(top10$gene, cc_markers))) 
  
  
  ###cluster 6 has the 'unclear' cells in this case 
  
  #removal of cluster of indistinct cells that overlap with 'unclear' ones 
  dim(pbmc)
  cells_to_remove <- WhichCells(pbmc, idents = c(6))
  pbmc <- subset(pbmc, cells = setdiff(Cells(pbmc), cells_to_remove))
  dim(pbmc)
  
  
  #repeat sctransform using subset object
  pbmc= SCTransform(pbmc) 
  dim(pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000, verbose = F)
  dim(pbmc)
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
  dim(pbmc)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.3,  )
  dim(pbmc)
  
  
  
  FeaturePlot(pbmc, features = cc_markers)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  
  
  
  #here and in 3051 i am adding a high rp + lowish s cluster to s but it could be left out 
  #nacl.t30.3004.cc = c("g1/S", "Stress", "M/G1", "G2/M", "S", "G2/M", "G2/M", "S", "G2/M", "M/G1", "G1")
  nacl.t30.3004.cc = c("G1_S", "Stress", "M_G1 + mating", "G2_M + mating", "S", "G2_M w/o mating", "G2_M + mating", "S", "Stress", "M_G1 w/o mating", "G1")
  
  
  
  
  names(nacl.t30.3004.cc) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, nacl.t30.3004.cc)
  
  #######################
  #find markers and make another heatmap using final (cc) clusters
  
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() -> top10
  
  DoHeatmap(pbmc, features = c(c(top10$gene, cc_markers), "MFA2", "CTT1")) 
  
  
  
  
  
  #############################################
  #########################
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
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
  dim(pbmc)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5,  )
  dim(pbmc)
  
  
  FeaturePlot(pbmc, features = cc_markers)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10, dims = c(2,3))
  
  
  #mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 15) %>%
    ungroup() -> top10
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers, "SPG1")) 
  
  #################################################
  
  
  
  #########################
  pbmc_data <- Read10X(data.dir = "/Users/noahalexander/SP_3051_rep1/SP_3051_rep1_t10//filtered_feature_bc_matrix/seurat_in/")
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
  pbmc= RunPCA(pbmc, npcs = 12, features = df$SGD)
  dim(pbmc)
  pbmc= FindNeighbors(pbmc,dims = 1:12, reduction = "pca") 
  dim(pbmc)
  pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5,  )
  dim(pbmc)
  
  
  FeaturePlot(pbmc, features = cc_markers)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10)
  DimPlot(pbmc, label = TRUE, pt.size = 1, label.size = 10, dims = c(2,3))
  
  
  #mat = GetAssayData(object = pbmc, assay = "SCT", slot = "counts")
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
  pbmc.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 24) %>%
    ungroup() -> top10
  
  DoHeatmap(pbmc, features = c(top10$gene, cc_markers, "SPG1"))
  
  
