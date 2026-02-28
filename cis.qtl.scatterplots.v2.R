pA <- ggplot(dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
geom_point(size = 1.2, alpha = 0.6, stroke = 0) +
scale_color_manual(
name = "Significance\n(FDR \u2264 0.05)",
values = c(
"Both timepoints" = "#D55E00",
"t30 only" = "#0072B2",
"t0 only" = "#009E73",
"Neither" = "grey70"
)
) +
labs(
x = "log(LOD) at t0 time point",
y = "log(LOD) at t30 time point"
) +
annotate("text",
x = -Inf, y = Inf,
label = "A",
hjust = -0.5, vjust = 1.5,
size = 6, fontface = "bold") +
theme_classic(base_size = 12) +
theme(
axis.title = element_text(size = 16),
axis.text = element_text(size = 11, color = "black"),
legend.title = element_text(size = 13),
legend.text = element_text(size = 12),
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

Then do the same for the stationary phase plot, but label it B and save as pB:

pB <- ggplot(dl, aes(log(LOD.x), log(LOD.y), color = vec)) +
geom_point(size = 1.2, alpha = 0.6, stroke = 0) +
scale_color_manual(
name = "Significance\n(FDR \u2264 0.05)",
values = c(
"Both timepoints" = "#D55E00",
"t10 only" = "#0072B2",
"t0 only" = "#009E73",
"Neither" = "grey70"
)
) +
labs(
x = "log(LOD) at t0 time point",
y = "log(LOD) at t10 time point"
) +
annotate("text",
x = -Inf, y = Inf,
label = "B",
hjust = -0.5, vjust = 1.5,
size = 6, fontface = "bold") +
theme_classic(base_size = 12) +
theme(
axis.title = element_text(size = 16),
axis.text = element_text(size = 11, color = "black"),
legend.title = element_text(size = 13),
legend.text = element_text(size = 12),
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

Now arrange them side by side using grid:

library(grid)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

print(pA, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pB, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
