#code to generate a histogram representing the number of contexts in which hotspots are observed 

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
BYxRM_NaCl_t0_G1 = unique(df$bin)
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
BYxRM_NaCl_t30_G1 = unique(df$bin)
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
BYxRM_NaCl_t30_M_G1_no_mating = unique(df$bin)
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
  BYxRM_NaCl_t0_combined=BYxRM_NaCl_t0_combined, BYxRM_NaCl_t0_G1=BYxRM_NaCl_t0_G1, BYxRM_NaCl_t0_G1_S=BYxRM_NaCl_t0_G1_S, BYxRM_NaCl_t0_G2_M_mating=BYxRM_NaCl_t0_G2_M_mating,
  BYxRM_NaCl_t0_M_G1_mating=BYxRM_NaCl_t0_M_G1_mating, BYxRM_NaCl_t30_combined=BYxRM_NaCl_t30_combined, 
  BYxRM_NaCl_t30_G1=BYxRM_NaCl_t30_G1,BYxRM_NaCl_t30_G1_S=BYxRM_NaCl_t30_G1_S, BYxRM_NaCl_t30_G2_M_mating=BYxRM_NaCl_t30_G2_M_mating,
  BYxRM_NaCl_t30_G2_M_no_mating=BYxRM_NaCl_t30_G2_M_no_mating,BYxRM_NaCl_t30_M_G1_no_mating=BYxRM_NaCl_t30_M_G1_no_mating, BYxRM_sp_t0_I=BYxRM_sp_t0_I, 
  BYxRM_sp_t0_II=BYxRM_sp_t0_II, BYxRM_sp_t0_IV=BYxRM_sp_t0_IV, BYxRM_sp_t10_I=BYxRM_sp_t10_I,BYxRM_sp_t10_II=BYxRM_sp_t10_II,
  BYxRM_sp_t10_III=BYxRM_sp_t10_III,BYxRM_sp_t10_IV=BYxRM_sp_t10_IV, BYxRM_sp_t10_combined=BYxRM_sp_t10_combined, BYxRM_sp_t0_combined=BYxRM_sp_t0_combined, 
  BYxRM_NaCl_t0_G2_M_no_mating=BYxRM_NaCl_t0_G2_M_no_mating, BYxRM_NaCl_t0_S=BYxRM_NaCl_t0_S, BYxRM_NaCl_t30_S = BYxRM_NaCl_t30_S ,  CBSxYJM_NaCl_t0_S = CBSxYJM_NaCl_t0_S)






unique_per_list <- lapply(test.input, unique)

all_elements <- unique(unlist(unique_per_list))

presence_counts <- sapply(all_elements, function(x) {
  sum(sapply(unique_per_list, function(vec) x %in% vec))
})

grouped_counts <- ifelse(presence_counts >= 5, "≥5", as.character(presence_counts))

levels_order <- c(as.character(1:4), "≥5")
grouped_factor <- factor(grouped_counts, levels = levels_order)

freq_table <- table(grouped_factor)

bp <- barplot(freq_table,
  xlab = "Number of State Contexts",
  ylab = "Number of Hotspots",
  col = "skyblue",
  cex.lab = 1.6,
  cex.axis = 1.2,
  ylim = c(0, max(freq_table) * 1.1)
)

text(x = bp, y = freq_table, labels = freq_table, pos = 3)


#######################################showing all state context numbers vs compressing above 5

levels_order <- as.character(1:max(presence_counts))
grouped_factor <- factor(as.character(presence_counts), levels = levels_order)

freq_table <- table(grouped_factor)

bp <- barplot(freq_table,
  xlab = "Number of State Contexts",
  ylab = "Number of Hotspots",
  col = "skyblue",
  cex.lab = 1.5,
  cex.axis = 1.3,
  ylim = c(0, max(freq_table) * 1.1)
)

nz <- freq_table > 0
text(x = bp[nz], y = freq_table[nz], labels = freq_table[nz], pos = 3)




########################################or using ggplot2
df_plot <- data.frame(
  n_contexts = factor(names(freq_table), levels = names(freq_table)),
  n_hotspots = as.integer(freq_table)
)

p <- ggplot(df_plot, aes(x = n_contexts, y = n_hotspots)) +
  geom_col(fill = "skyblue", color = "black", linewidth = 0.2) +
  geom_text(aes(label = n_hotspots), vjust = -0.3, size = 5) +
  labs(
    x = "Number of State Contexts",
    y = "Number of Hotspots"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(0, max(df_plot$n_hotspots) * 1.1)) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    plot.margin = margin(10, 10, 10, 20)
  )

print(p)

ggsave("hotspot_context_hist.pdf", p, width = 7, height = 4.5, useDingbats = FALSE)
