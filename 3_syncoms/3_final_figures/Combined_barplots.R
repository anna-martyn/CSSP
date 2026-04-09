# Set working directory and load packages --------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pkg <- c("data.table", "ggplot2", "ggtext", "ggh4x")
for(pk in pkg){
  library(pk, character.only = T)
}

# Settings ---------------------------------------------------------------------
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size = 6, color="black"),
  legend.text=element_text(size = 6, color="black"),
  legend.key=element_blank(),
  axis.title.y=element_text(size = 6),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size = 6, color="black")
)

colors_geno <- c(
  "WT"     = "#A9C289",
  "symrk"  = "#FEDA8B",
  "ccamk"  = "#FDB366",
  "nsp1"   = "#C0E4EF",
  "nsp2"   = "#6EA6CD"
)

genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Load data --------------------------------------------------------------------
df_ord_Lj <- fread(
  input = "../1_Lotus/7_stackedbp_barplots/3_tables/df_ord_Lj.csv",
  drop = 1
)
df_ord_Hv <- fread(
  input = "../2_Hordeum/6_stackedbp_barplots/3_tables/df_ord_Hv.csv",
  drop = 1
)
df_ord <- rbind(df_ord_Lj, df_ord_Hv)
df_ord[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Genotype = factor(Genotype, levels = names(colors_geno))
)]

# Visualisation ----------------------------------------------------------------
p_sig <- ggplot(df_ord, aes(x = Order, y = mean_RA, fill = Genotype)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    width = 0.8,
    alpha = 0.9
  ) +
  geom_errorbar(
    aes(ymin = mean_RA - 1.96 * se, ymax = mean_RA + 1.96 * se),
    width = 0.3,
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    data = df_ord,
    aes(
      x = Order,
      y = mean_RA + 1.96 * se + 0.015,
      label = Sig,
      fill = Genotype
    ),
    position = position_dodge(width = 0.9),
    inherit.aes = FALSE,
    size = 4
  ) +
  # facet_wrap(~Host + Compartment, scales = "free_x", nrow = 1) +
  facet_wrap2(
    vars(Host, Compartment),
    strip = strip_nested(),
    nrow = 1,
    scales = "free_x"
  ) +
  scale_fill_manual(values = colors_geno, labels = genotype_labels_legend) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "", y = "Relative Abundance", fill = "Genotype") +
  main_theme +
  force_panelsizes(cols = c(5, 4, 8, 5), rows = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_markdown(),
    strip.text = element_text(face = "bold", size = rel(1)),
    legend.position = "none",
    plot.title = element_text(size = 6)
  )

ggsave(
  filename = "2_temp_figures/Plot_RA_sign_orders.pdf",
  plot = p_sig,
  width = 21,
  height = 6,
  unit = "cm"
)
saveRDS(object = p_sig, file = "1_rds_files/Plot_RA_sign_orders.rds")
