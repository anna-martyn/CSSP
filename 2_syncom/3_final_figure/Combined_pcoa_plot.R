# Set working directory and load packages --------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(ggh4x)
library(ggtext)
library(multcompView)
library(data.table)

# Settings ---------------------------------------------------------------------
# Set colors for genotypes.
colors <- c("WT" = "#A9C289", "symrk" = "#FEDA8B", "ccamk" = "#FDB366",
            "nsp1" = "#C0E4EF", "nsp2" = "#6EA6CD")

legend_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Set main theme for plot.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text.x = element_text(size = 8, colour = "black"),
  axis.text.y = element_text(size = 8, colour = "black"),
  legend.background=element_blank(),
  legend.key=element_blank(),
  text=element_text(family="sans")
)

# Load data --------------------------------------------------------------------
points_Lj <- fread("../1_Lotus/6_cpcoa_pcoa/pcoa/PCoA_points_Lj.csv")
points_Hv <- fread("../2_Hordeum/5_cpcoa_pcoa/pcoa/PCoA_points_Hv.csv")

segments_Lj <- fread("../1_Lotus/6_cpcoa_pcoa/pcoa/PCoA_segments_Lj.csv")
segments_Hv <- fread("../2_Hordeum/5_cpcoa_pcoa/pcoa/PCoA_segments_Hv.csv")

text_dt_Lj <- fread("../1_Lotus/6_cpcoa_pcoa/pcoa/PCoA_text_Lj.csv")
text_dt_Hv <- fread("../2_Hordeum/5_cpcoa_pcoa/pcoa/PCoA_text_Hv.csv")

points <- rbind(points_Lj[,-8], points_Hv[,-6])
segments <- rbind(segments_Lj[,-8], segments_Hv[,-6])
text_dt <- rbind(text_dt_Lj, text_dt_Hv)

points[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Genotype = factor(Genotype, levels = names(colors))
)]

segments[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Genotype = factor(Genotype, levels = names(colors))
)]

text_dt[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  text = gsub("-", "\n", text)
)]

text_dt[,text:=gsub("52", "52.0", text)]

# Visualisation ----------------------------------------------------------------
p <- ggplot(points, aes(x=x, y=y, color=Genotype)) +
  geom_point(size=1.5, alpha=0.7) +
  facet_grid(Host~Compartment, switch = "y")+
  geom_segment(data=segments, aes(x = x, y = y, xend = seg_x,
                                  yend = seg_y, color = Genotype),
               alpha=0.5) +
  geom_label(data = text_dt, aes(x = -0.175, y = -0.25, label = text),
             colour = "black", fill = "grey", alpha = 0.2, size = 8/.pt)+
  scale_color_manual(values=colors, labels=legend_labels) +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  labs(
    x = "PCo 1",
    y = "PCo 2"
  ) +
  main_theme +
  theme(
    plot.title = element_text(face="bold", size=8, hjust=0),
    legend.text = element_markdown(size=8, color="black"),
    strip.text = element_text(size = 8, colour = "black", face = "bold"),
    strip.placement = "outside",
    axis.title.x=element_text(size = 8, colour = "black"),
    axis.title.y=element_text(size = 8, colour = "black"),
    legend.key.size = unit(0.25, "cm")
  )+
  NULL
p

ggsave("PCoA_plot_matched.pdf", p, width = 5, height = 13, unit = "cm")
saveRDS(p, file = "PCoA_plot_matched.rds")
