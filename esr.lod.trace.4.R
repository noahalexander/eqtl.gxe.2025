library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

# ---------- Load + reshape helper ----------
prep_lod_long <- function(path, time_label) {
    lod <- readRDS(path)$LOD
    lod <- as.data.frame(t(lod))
    lod$marker <- rownames(lod)
    lod$pos <- seq_len(nrow(lod))
    
    colnames(lod)[1:3] <- c("RP", "iESR", "RiBi")
    
    lod %>%
        select(marker, pos, RP, iESR, RiBi) %>%
        pivot_longer(cols = c(RP, iESR, RiBi),
            names_to = "state",
            values_to = "LOD") %>%
        mutate(time = time_label)
}

# ---------- Load BYxRM ----------
A_t0 <- prep_lod_long(
    "/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/state_pheno_LODs.RDS",
    "t0"
)

A_t30 <- prep_lod_long(
    "/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/state_pheno_LODs.RDS",
    "t30"
)

plot_df <- bind_rows(A_t0, A_t30) %>%
    mutate(
        state = factor(state, levels = c("RP", "iESR", "RiBi")),
        time  = factor(time, levels = c("t0", "t30"))
    )

# ---------- Chromosome boundaries (your method) ----------
chr_map <- A_t30 %>%
    distinct(marker, pos) %>%
    mutate(chr = str_extract(marker, "^[^_]+")) %>%
    arrange(pos)

rle_chr <- rle(chr_map$chr)
cuts <- c(0, cumsum(rle_chr$lengths))
breaks_all <- cuts[-length(cuts)]   # includes 0 (needed so breaks/labels lengths match)
labels_all <- rle_chr$values
boundary_lines <- breaks_all[breaks_all > 0]

# ---------- Softer palette (what you’re currently using) ----------
state_colors <- c(
    "RP"   = "#6BA292",  # soft sage
    "iESR" = "#C77C8A",  # dusty rose
    "RiBi" = "#6C7A89"   # steel gray-blue
)

# ---------- Y cap ----------
y_cap <- 75

# ---------- Plot ----------
p <- ggplot(plot_df, aes(x = pos, y = LOD, color = state, linetype = time)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 4, linewidth = 0.5, color = "grey40") +
    geom_vline(xintercept = boundary_lines, linetype = "dashed",
        linewidth = 0.35, color = "grey92") +
    
    # keep the cap line (as a "cap" indicator), but not as a border
    geom_hline(yintercept = y_cap, linewidth = 0.4, color = "grey25") +
    
    scale_color_manual(values = state_colors, name = NULL) +
    scale_linetype_manual(values = c("t0" = "solid", "t30" = "dashed"), name = NULL) +
    scale_x_continuous(breaks = breaks_all, labels = labels_all, expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_cap), expand = c(0, 0), oob = squish) +
    labs(x = NULL, y = "LOD") +
    theme_classic(base_size = 16) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        
        # make all four sides identical (instead of mixing geoms + axis lines)
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.6),
        
        # keep legend unboxed
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        
        plot.margin = margin(5, 10, 5, 5)
    )

p


ggsave(
    filename = "Figure_esr_QTL.png",
    plot = p,
    width = 12,
    height = 6,
    dpi = 400
)
