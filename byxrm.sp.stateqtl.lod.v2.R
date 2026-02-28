# ---------------- SP t0 ----------------
t0 <- prep_cc_long(
    path_cc_lod = "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/cell_cycle_assignment_LOD.RDS",
    path_state_pheno_lods = "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/state_pheno_LODs.RDS",
    time_label = "t0"
)

# colors (match your base-R intent; tweak as desired)
t0_levels <- unique(t0$data$series)
t0_colors <- c(
    "t0 I"   = "lightgreen",
    "t0 II"  = "darkgreen",
    "t0 III" = "lightblue",
    "t0 IV"  = "darkblue",
    "t0 V"   = "black",
    "t0 VI"  = "darkred"
)

p_t0 <- ggplot(t0$data, aes(x = pos, y = LOD, color = series)) +
    geom_line(size = 0.7) +
    geom_hline(yintercept = 4, size = 0.6, color = "grey40") +
    geom_vline(xintercept = t0$breaks, linetype = "dashed", size = 0.4, color = "grey75") +
    scale_color_manual(values = t0_colors, breaks = names(t0_colors), name = NULL) +
    scale_x_continuous(breaks = t0$breaks, labels = t0$chr_labels, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 50)) +
    labs(x = NULL, y = "LOD") +
    pub_theme()

p_t0

ggsave(
    filename = "Figure4.1_QTL.png",
    plot = p_t0,
    width = 12,
    height = 6,
    dpi = 400
)


# ---------------- SP t10 ----------------
t10 <- prep_cc_long(
    path_cc_lod = "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/cell_cycle_assignment_LOD.RDS",
    path_state_pheno_lods = "/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/state_pheno_LODs.RDS",
    time_label = "t10"
)

t10_colors <- c(
    "t10 I"   = "lightgreen",
    "t10 II"  = "darkgreen",
    "t10 III" = "lightblue",
    "t10 IV"  = "darkblue",
    "t10 V"   = "black",
    "t10 VI"  = "darkred",
    "t10 VII" = "purple"
)

p_t10 <- ggplot(t10$data, aes(x = pos, y = LOD, color = series)) +
    geom_line(size = 0.7) +
    geom_hline(yintercept = 4, size = 0.6, color = "grey40") +
    geom_vline(xintercept = t10$breaks, linetype = "dashed", size = 0.4, color = "grey75") +
    scale_color_manual(values = t10_colors, breaks = names(t10_colors), name = NULL) +
    scale_x_continuous(breaks = t10$breaks, labels = t10$chr_labels, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 50)) +
    labs(x = NULL, y = "LOD") +
    pub_theme()

p_t10
ggsave(
    filename = "Figure4.2_QTL.png",
    plot = p_t10,
    width = 12,
    height = 6,
    dpi = 400
)


########

combined_plot <- p_t0 / p_t10

combined_plot

ggsave(
    filename = "Figure6_SP_state_QTL.png",
    plot = combined_plot,
    width = 12,
    height = 6,
    dpi = 400
)
