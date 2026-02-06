

library(dplyr)
library(ggplot2)

df.nt0 = readRDS("repeat_fine_mapping/combined/A/nacl/t0/cis_only_test_CombinedResults.RDS")
df.nt0 = df.nt0$combined
#df = subset(df, df$FDR <= 0.05)
dim(df.nt0)

df.3 = readRDS("repeat_fine_mapping/combined/A/nacl/t30/cis_only_test_CombinedResults.RDS")
df.3 = df.3$combined
#df.3 = subset(df.3, df.3$FDR <= 0.05)
dim(df.3)


#dd = merge(df, df.3, by = "transcript", all=T)


df.spt0 = readRDS("repeat_fine_mapping/combined/A/sp/t0/cis_only_test_CombinedResults.RDS")
df.spt0= df.spt0$combined
#df.spt0 = subset(df.spt0, df.spt0$FDR <= 0.05)



df.t10 = readRDS("repeat_fine_mapping/combined/A/sp/t10/cis_only_test_CombinedResults.RDS")
df.t10 = df.t10$combined
#df.t10 = subset(df.t10, df.t10$FDR <= 0.05)

dl = merge(df.spt0, df.t10, by = "transcript", all=T)


dl <- merge(df.nt0, df.3, by = "transcript", all = TRUE)



############BYxRM salt perturbation 

dl = merge(df.nt0, df.3, by = "transcript", all=T)

dl$FDR.x[is.na(dl$FDR.x)] <- 1
dl$FDR.y[is.na(dl$FDR.y)] <- 1

vec= vector()
for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <= 0.05 & dl$FDR.y[i] <= 0.05) {
    vec[i] = "both"
  }
}
for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] >0.05 & dl$FDR.y[i] <= 0.05) {
    vec[i] = "t30.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <=0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "t0.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] > 0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "neither"
  }
}

dl$vec = vec

dl$vec <- factor(
  dl$vec,
  levels = c("both", "t30.only", "t0.only", "neither"),
  labels = c(
    "Both timepoints",
    "t30 only",
    "t0 only",
    "Neither"
  )
)




p <- ggplot(dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point(size = 1.2, alpha = 0.6, stroke = 0) +
  scale_color_manual(
    name = "Significance\n(FDR \u2264 0.05)",
    values = c(
      "Both timepoints" = "#D55E00",  # orange-red
      "t30 only"        = "#0072B2",  # blue
      "t0 only"         = "#009E73",  # green
      "Neither"         = "grey70"
    )
  ) +
  labs(
    x = "log(LOD) at t0 time point",
    y = "log(LOD) at t30 time point"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    legend.key = element_blank(),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(
      fill = "white",
      color = "grey85",
      size = 0.3
    ),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

p


############# BYxRM stationary phase



dl = merge(df.spt0, df.t10, by = "transcript", all=T)

dl$FDR.x[is.na(dl$FDR.x)] <- 1
dl$FDR.y[is.na(dl$FDR.y)] <- 1

vec= vector()
for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <= 0.05 & dl$FDR.y[i] <= 0.05) {
    vec[i] = "Both timepoints"
  }
}
for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] >0.05 & dl$FDR.y[i] <= 0.05) {
    vec[i] = "t10 only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <=0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "t0 only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] > 0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "Neither"
  }
}

dl$vec = factor(vec)

dl$vec <- factor(
  dl$vec,
  levels = c("Both timepoints", "t10 only", "t0 only", "Neither"),
  labels = c(
    "Both timepoints",
    "t10 only",
    "t0 only",
    "Neither"
  )
)



p <- ggplot(dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point(size = 1.2, alpha = 0.6, stroke = 0) +
  scale_color_manual(
    name = "Significance\n(FDR \u2264 0.05)",
    values = c(
      "Both timepoints" = "#D55E00",  # orange-red
      "t10 only"        = "#0072B2",  # blue
      "t0 only"         = "#009E73",  # green
      "Neither"         = "grey70"
    )
  ) +
  labs(
    x = "log(LOD) at t0 time point",
    y = "log(LOD) at t10 time point"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    legend.key = element_blank(),
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(
      fill = "white",
      color = "grey85",
      size = 0.3
    ),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

p
