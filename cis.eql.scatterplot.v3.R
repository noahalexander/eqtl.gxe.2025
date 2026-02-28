# ============================================================
# Two-panel scatterplots (A,B) with ONE shared legend on right
# Legend title bold + size 15; legend has no box
# Requires: df.nt0, df.3, df.spt0, df.t10 already loaded as in your script
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ----------------------------
# Panel A data: NaCl t0 vs t30
# ----------------------------
dl_nt <- merge(df.nt0, df.3, by = "transcript", all = TRUE)
dl_nt$FDR.x[is.na(dl_nt$FDR.x)] <- 1
dl_nt$FDR.y[is.na(dl_nt$FDR.y)] <- 1

dl_nt$vec <- case_when(
  dl_nt$FDR.x <= 0.05 & dl_nt$FDR.y <= 0.05 ~ "Both timepoints",
  dl_nt$FDR.x >  0.05 & dl_nt$FDR.y <= 0.05 ~ "Perturbation only",
  dl_nt$FDR.x <= 0.05 & dl_nt$FDR.y >  0.05 ~ "Baseline only",
  TRUE                                      ~ "Neither"
)

# ----------------------------
# Panel B data: SP t0 vs t10
# ----------------------------
dl_sp <- merge(df.spt0, df.t10, by = "transcript", all = TRUE)
dl_sp$FDR.x[is.na(dl_sp$FDR.x)] <- 1
dl_sp$FDR.y[is.na(dl_sp$FDR.y)] <- 1

dl_sp$vec <- case_when(
  dl_sp$FDR.x <= 0.05 & dl_sp$FDR.y <= 0.05 ~ "Both timepoints",
  dl_sp$FDR.x >  0.05 & dl_sp$FDR.y <= 0.05 ~ "Perturbation only",
  dl_sp$FDR.x <= 0.05 & dl_sp$FDR.y >  0.05 ~ "Baseline only",
  TRUE                                      ~ "Neither"
)

# ----------------------------
# Common factor levels + scale
# ----------------------------
common_levels <- c("Both timepoints", "Perturbation only", "Baseline only", "Neither")
dl_nt$vec <- factor(dl_nt$vec, levels = common_levels)
dl_sp$vec <- factor(dl_sp$vec, levels = common_levels)

sig_scale <- scale_color_manual(
  name   = "Significance\n(FDR ≤ 0.05)",
  breaks = common_levels,
  limits = common_levels,
  drop   = FALSE,
  values = c(
    "Both timepoints"   = "#D55E00",
    "Perturbation only" = "#0072B2",
    "Baseline only"     = "#009E73",
    "Neither"           = "grey70"
  )
)

# ----------------------------
# Build Panel A plot (no legend)
# ----------------------------
pA <- ggplot(dl_nt, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point(size = 1.2, alpha = 0.6) +
  sig_scale +
  labs(
    x = "log(LOD) at t0 time point",
    y = "log(LOD) at t30 time point"
  ) +
  annotate("text", x = -Inf, y = Inf, label = "A",
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold") +
  theme_classic(base_size = 12) +
  theme(
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

# ----------------------------
# Build Panel B plot (legend kept)
# ----------------------------
pB <- ggplot(dl_sp, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point(size = 1.2, alpha = 0.6) +
  sig_scale +
  labs(
    x = "log(LOD) at t0 time point",
    y = "log(LOD) at t10 time point"
  ) +
  annotate("text", x = -Inf, y = Inf, label = "B",
           hjust = -0.5, vjust = 1.5, size = 6, fontface = "bold") +
  theme_classic(base_size = 12) +
  theme(
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 11, color = "black"),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

# ----------------------------
# Combine with ONE shared legend on the right
# and style legend title to match figure
# ----------------------------
final_plot <- (pA + pB) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 12)
  )

final_plot

# Optional export (PDF recommended for papers)
ggsave("figure_AB_shared_legend.png", final_plot, width = 12, height = 6, units = "in")
