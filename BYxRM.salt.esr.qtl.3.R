library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

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
    mutate(time = time_label,
           series = paste(time, state))
}

# ---------- Load data ----------
df_t0  <- prep_lod_long(
  "/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/state_pheno_LODs.RDS",
  "t0"
)

df_t30 <- prep_lod_long(
  "/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/state_pheno_LODs.RDS",
  "t30"
)

plot_df <- bind_rows(df_t0, df_t30)

# ---------- Chromosome boundaries ----------
chr_map <- df_t30 %>%
  distinct(marker, pos) %>%
  mutate(chr = str_extract(marker, "^[^_]+")) %>%
  arrange(pos)

rle_chr <- rle(chr_map$chr)
cuts <- c(0, cumsum(rle_chr$lengths))
breaks <- cuts[-length(cuts)]
labels <- rle_chr$values

# ---------- Colors ----------
series_colors <- c(
  "t0 RP"    = "lightgreen",
  "t0 iESR"  = "darkgreen",
  "t0 RiBi"  = "lightblue",
  "t30 RP"   = "orange",
  "t30 iESR" = "red",
  "t30 RiBi" = "yellow"
)

# ---------- Plot ----------
p <- ggplot(plot_df, aes(x = pos, y = LOD, color = series)) +
  geom_line(size = 0.6) +
  geom_hline(yintercept = 4, size = 0.5, color = "grey40") +
  geom_vline(xintercept = breaks, linetype = "dashed",
             size = 0.4, color = "grey75") +
  scale_color_manual(values = series_colors,
                     breaks = names(series_colors),
                     name = NULL) +
  scale_x_continuous(breaks = breaks,
                     labels = labels,
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 75)) +
  labs(x = NULL, y = "LOD") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),

    # Full box around plotting panel
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),

    # Legend outside with box
    legend.position = "right",
    legend.background = element_rect(color = "black", fill = "white", size = 0.6),
    legend.box.background = element_rect(color = "black", fill = NA, size = 0.6),

    legend.key = element_blank(),
    legend.text = element_text(size = 14),
    axis.title.y = element_text(size = 18),

    plot.margin = margin(10, 20, 10, 10)
  )

p
