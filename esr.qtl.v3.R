library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(grid)

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

# ---------- YOUR chr boundary method (breaks+labels lengths match) ----------
get_chr_meta <- function(df_t30) {
  chr_map <- df_t30 %>%
    distinct(marker, pos) %>%
    mutate(chr = str_extract(marker, "^[^_]+")) %>%
    arrange(pos)

  rle_chr <- rle(chr_map$chr)
  cuts <- c(0, cumsum(rle_chr$lengths))

  breaks_all <- cuts[-length(cuts)]   # includes 0
  labels_all <- rle_chr$values        # same length

  list(
    vlines  = breaks_all[breaks_all > 0],
    xbreaks = breaks_all,
    xlabels = labels_all
  )
}

# ---------- Aesthetics ----------
state_colors <- c(
  "RP"   = "#2E7D32",
  "iESR" = "#D32F2F",   # red
  "RiBi" = "#111111"
)
y_cap <- 75

make_plot <- function(plot_df, chr_meta, panel_letter, show_x = TRUE) {
  g <- ggplot(plot_df, aes(x = pos, y = LOD, color = state, linetype = time)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 4, linewidth = 0.5, color = "grey40") +
    geom_vline(xintercept = chr_meta$vlines, linetype = "dashed",
               linewidth = 0.35, color = "grey85") +
    geom_hline(yintercept = y_cap, linewidth = 0.4, color = "grey60") +
    scale_color_manual(values = state_colors, name = NULL) +
    scale_linetype_manual(values = c("t0" = "solid", "t30" = "dashed"), name = NULL) +
    scale_x_continuous(breaks = chr_meta$xbreaks,
                       labels = chr_meta$xlabels,
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_cap), expand = c(0, 0), oob = squish) +
    labs(x = NULL, y = paste0("LOD (capped at ", y_cap, ")")) +
    annotate("text",
             x = -Inf, y = Inf, label = panel_letter,
             hjust = -0.4, vjust = 1.4,
             fontface = "bold", size = 6) +
    theme_classic(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.border = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.position = "right",
      legend.key = element_blank(),
      plot.margin = margin(5, 10, 5, 22)
    )

  if (!show_x) {
    g <- g + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }
  g
}

# =======================
# BYxRM (A)
# =======================
A_t0  <- prep_lod_long("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/state_pheno_LODs.RDS",  "t0")
A_t30 <- prep_lod_long("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/state_pheno_LODs.RDS", "t30")

plot_A <- bind_rows(A_t0, A_t30) %>%
  mutate(
    state = factor(state, levels = c("RP", "iESR", "RiBi")),
    time  = factor(time,  levels = c("t0", "t30"))
  )

meta_A <- get_chr_meta(A_t30)
pA <- make_plot(plot_A, meta_A, "A", show_x = FALSE)

# =======================
# CBSxYJM (3004) (B)
# =======================
B_t0  <- prep_lod_long("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/state_pheno_LODs.RDS",  "t0")
B_t30 <- prep_lod_long("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/state_pheno_LODs.RDS", "t30")

plot_B <- bind_rows(B_t0, B_t30) %>%
  mutate(
    state = factor(state, levels = c("RP", "iESR", "RiBi")),
    time  = factor(time,  levels = c("t0", "t30"))
  )

meta_B <- get_chr_meta(B_t30)
pB <- make_plot(plot_B, meta_B, "B", show_x = TRUE)

# =======================
# DRAW BOTH IN THE SAME PLOTTING WINDOW (NO SAVING)
# =======================
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))

print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
