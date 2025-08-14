library(ggplot2)
df.0 = readRDS("repeat_fine_mapping/combined/3004/nacl/t0/cis_only_test_CombinedResults.RDS")
df.0= df.0$combined
#df.spt0 = subset(df.spt0, df.spt0$FDR <= 0.05)

df.30 = readRDS("repeat_fine_mapping/combined/3004//nacl/t30/cis_only_test_CombinedResults.RDS")
df.30 = df.30$combined

dl = merge(df.0, df.30, by = "transcript", all=T)

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
    vec[i] = "nacl.t0.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] <=0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "nacl.t30.only"
  }}

for (i in 1:nrow(dl)) {
  if (dl$FDR.x[i] > 0.05 & dl$FDR.y[i] > 0.05) {
    vec[i] = "neither"
  }
}

dl$vec = vec
ggplot(data = dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point() +
  scale_color_manual(values = c("both" = "red", "nacl.t0.only" = "blue", "nacl.t30.only" = "green", "neither" = "black"))

ggplot(data = dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
  geom_point() +
  scale_color_manual(values = c("both" = "red", "nacl.t0.only" = "blue", "nacl.t30.only" = "green", "neither" = "black")) +
  labs(x = "log(LOD)", y = "log(LOD)" )+
  theme(
    axis.title = element_text(size = 18),      # Axis label size
    axis.text = element_text(size = 14),       # Axis tick label size
    legend.text = element_text(size = 16),     # Legend text size
    legend.title = element_text(size = 18)     # Legend title size
  )

