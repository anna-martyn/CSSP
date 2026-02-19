# Set working directory and load packages --------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(data.table)
library(ggtext)
library(ggh4x)

# Settings ---------------------------------------------------------------------
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size=8, color="black"),
  legend.text=element_text(size=8, color="black"),
  legend.key=element_blank(),
  axis.title.y=element_text(size=8),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size=8, color="black")
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

dodge <- position_dodge(width = 0.9)

# Load data --------------------------------------------------------------------
df_ord_Lj <- fread("../1_Lotus/7_stackedbp_barplots/df_ord_Lj.csv", drop = 1)
df_ord_Hv <- fread("../2_Hordeum/6_stackedbp_barplots/df_ord_Hv.csv", drop = 1)
df_ord <- rbind(df_ord_Lj, df_ord_Hv)
df_ord[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Genotype = factor(Genotype, levels = names(colors_geno))
)]

# Visualisation ----------------------------------------------------------------
p_sig <- ggplot(df_ord, aes(x = Order, y = mean_RA, fill = Genotype)) +
  geom_bar(stat = "identity", position = dodge, width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin=mean_RA-1.96*se, ymax=mean_RA+1.96*se),
                width = 0.3, position = dodge) +
  geom_text(
    data = df_ord,
    aes(x = Order, y=mean_RA+1.96*se+0.015, label=Sig, fill=Genotype),
    position=dodge,
    inherit.aes=FALSE,
    size=4
  ) + 
  # facet_wrap(~Host + Compartment, scales = "free_x", nrow = 1) +
  facet_wrap2(vars(Host, Compartment), strip = strip_nested(), nrow = 1,
              scales = "free_x")+
  scale_fill_manual(values = colors_geno, labels=genotype_labels_legend) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x="", y="Relative Abundance", fill="Genotype") +
  main_theme +
  force_panelsizes(cols = c(5, 4, 8, 5), rows = 1)+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    legend.text = element_markdown(),
    strip.text = element_text(face="bold", size=rel(1)),
    legend.position = "none",
    plot.title = element_text(size = 8)
  )
p_sig

ggsave("Plot_RA_sign_orders.pdf", p_sig,
       width = 21, height = 6, unit = "cm")
saveRDS(p_sig, file="Plot_RA_sign_orders.rds")
