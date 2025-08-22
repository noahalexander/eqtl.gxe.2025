library(ComplexHeatmap)
library(ggplot2)
######3004 t0
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t0/hotspot_peaks.RDS")
length(df)
df = df$combined
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t0_combined = unique(df$bin)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t0/hotspot_peaks.RDS")
length(df)
df = df$G1_S
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t0_G1_S = unique(df$bin)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t0/hotspot_peaks.RDS")
length(df)
df = df$G2_M_mating
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t0_G2_M_mating = unique(df$bin)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t0/hotspot_peaks.RDS")
length(df)
df = df$M_G1_mating
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t0_M_G1_mating = unique(df$bin)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t0/hotspot_peaks.RDS")
length(df)
df = df$S
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t0_S = unique(df$bin)

######3004 t30
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t30/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t30_combined = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t30/hotspot_peaks.RDS")
df = df$G1_S
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t30_G1_S = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_mating
  df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t30_G2_M_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t30/hotspot_peaks.RDS")
df = df$M_G1_mating
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t30_M_G1_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t30/hotspot_peaks.RDS")
df = df$S
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t30_S = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004//nacl/t30/hotspot_peaks.RDS")
df = df$Stress
df = subset(df, df$in.hotspot == "TRUE")
CBSxYJM_NaCl_t30_Stress = unique(df$bin)


######3051 t0
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t0/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t0_combined = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t0/hotspot_peaks.RDS")
df = df$G1
df = subset(df, df$in.hotspot == "TRUE")
BYxRMNaCl_t0_G1 = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t0/hotspot_peaks.RDS")
df = df$G1_S
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t0_G1_S = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t0_G2_M_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_no_mating
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t0_G2_M_no_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t0/hotspot_peaks.RDS")
df = df$M_G1_mating
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t0_M_G1_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t0/hotspot_peaks.RDS")
df = df$S
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t0_S = unique(df$bin)

######3051 t30
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t30/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t30_combined = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t30/hotspot_peaks.RDS")
df = df$G1
df = subset(df, df$in.hotspot == "TRUE")
BYxRMNaCl_t30_G1 = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t30/hotspot_peaks.RDS")
df = df$G1_S
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t30_G1_S = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_mating
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t30_G2_M_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_no_mating
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t30_G2_M_no_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t30/hotspot_peaks.RDS")
df = df$M_G1_no_mating
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t30_M_G1_mating = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//nacl/t30/hotspot_peaks.RDS")
df = df$S
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_NaCl_t30_S = unique(df$bin)


######3051 s t0
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t0/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t0_combined = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t0/hotspot_peaks.RDS")
df = df$I
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t0_I = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t0/hotspot_peaks.RDS")
df = df$II
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t0_II = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t0/hotspot_peaks.RDS")
df = df$IV
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t0_IV = unique(df$bin)

######3051 s t10
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t10/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t10_combined = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t10/hotspot_peaks.RDS")
df = df$I
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t10_I = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t10/hotspot_peaks.RDS")
df = df$II
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t10_II = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t10/hotspot_peaks.RDS")
df = df$III
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t10_III = unique(df$bin)
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A//sp/t10/hotspot_peaks.RDS")
df = df$IV
df = subset(df, df$in.hotspot == "TRUE")
BYxRM_sp_t10_IV = unique(df$bin)


test.input = list(CBSxYJM_NaCl_t0_combined = CBSxYJM_NaCl_t0_combined, CBSxYJM_NaCl_t0_G1_S = CBSxYJM_NaCl_t0_G1_S, 
  CBSxYJM_NaCl_t0_G2_M_mating = CBSxYJM_NaCl_t0_G2_M_mating, CBSxYJM_NaCl_t0_M_G1_mating = CBSxYJM_NaCl_t0_M_G1_mating, 
  CBSxYJM_NaCl_t30_combined = CBSxYJM_NaCl_t30_combined, CBSxYJM_NaCl_t30_G1_S = CBSxYJM_NaCl_t30_G1_S, CBSxYJM_NaCl_t30_G2_M_mating=CBSxYJM_NaCl_t30_G2_M_mating,
  CBSxYJM_NaCl_t30_M_G1_mating=CBSxYJM_NaCl_t30_M_G1_mating, CBSxYJM_NaCl_t30_S=CBSxYJM_NaCl_t30_S, CBSxYJM_NaCl_t30_Stress=CBSxYJM_NaCl_t30_Stress, 
  BYxRM_NaCl_t0_combined=BYxRM_NaCl_t0_combined, BYxRMNaCl_t0_G1=BYxRMNaCl_t0_G1, BYxRM_NaCl_t0_G1_S=BYxRM_NaCl_t0_G1_S, BYxRM_NaCl_t0_G2_M_mating=BYxRM_NaCl_t0_G2_M_mating,
BYxRM_NaCl_t0_M_G1_mating=BYxRM_NaCl_t0_M_G1_mating, BYxRM_NaCl_t30_combined=BYxRM_NaCl_t30_combined, 
  BYxRMNaCl_t30_G1=BYxRMNaCl_t30_G1,BYxRM_NaCl_t30_G1_S=BYxRM_NaCl_t30_G1_S, BYxRM_NaCl_t30_G2_M_mating=BYxRM_NaCl_t30_G2_M_mating,
  BYxRM_NaCl_t30_G2_M_no_mating=BYxRM_NaCl_t30_G2_M_no_mating,BYxRM_NaCl_t30_M_G1_mating=BYxRM_NaCl_t30_M_G1_mating, BYxRM_sp_t0_I=BYxRM_sp_t0_I, 
  BYxRM_sp_t0_II=BYxRM_sp_t0_II, BYxRM_sp_t0_IV=BYxRM_sp_t0_IV, BYxRM_sp_t10_I=BYxRM_sp_t10_I,BYxRM_sp_t10_II=BYxRM_sp_t10_II,
  BYxRM_sp_t10_III=BYxRM_sp_t10_III,BYxRM_sp_t10_IV=BYxRM_sp_t10_IV, BYxRM_sp_t10_combined=BYxRM_sp_t10_combined, BYxRM_sp_t0_combined=BYxRM_sp_t0_combined)






unique_per_list <- lapply(test.input, unique)

all_elements <- unique(unlist(unique_per_list))

presence_counts <- sapply(all_elements, function(x) {
  sum(sapply(unique_per_list, function(vec) x %in% vec))
})

grouped_counts <- ifelse(presence_counts >= 5, "≥5", as.character(presence_counts))

levels_order <- c(as.character(1:4), "≥5")
grouped_factor <- factor(grouped_counts, levels = levels_order)

freq_table <- table(grouped_factor)

barplot(freq_table,
  xlab = "Number of State Contexts",
  ylab = "Number of Hotspots",
  col = "skyblue
