# Set working directory and load packages --------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(ggh4x)
library(multcompView)
library(data.table)

# Settings ---------------------------------------------------------------------
# Set colors for genotypes.
colors <- c("WT" = "#A9C289", "symrk" = "#FEDA8B", "ccamk" = "#FDB366",
            "nsp1" = "#C0E4EF", "nsp2" = "#6EA6CD")

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
weight_Lj <- fread("../1_Lotus/2_shoot_fw/LotusSC_shoot_fw_ANOVA.csv", drop = 1)
weight_Hv <- fread("../2_Hordeum/2_shoot_fw/HordeumSC_shoot_fw_ANOVA.csv", drop = 1)

weight_summary_Lj <- fread("../1_Lotus/2_shoot_fw/LotusSC_shoot_fw_significance_letters.csv",
                           drop = 1)
weight_summary_Hv <- fread("../2_Hordeum/2_shoot_fw/HordeumSC_shoot_fw_significance_letters.csv",
                           drop = 1)

weight_summary <- rbind(weight_summary_Lj, weight_summary_Hv)
weight <- rbind(weight_Lj[,-3], weight_Hv[,-(3:4)])

weight[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Genotype = factor(Genotype, levels = names(colors))
)]

weight_summary[,":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Genotype = factor(Genotype, levels = names(colors))
)]


# Visualisation ----------------------------------------------------------------
weight_summary$y_pos[c(6, 7, 10)] <- c(0.4, 0.38, 0.3)
weight_summary$y_pos[3:5] <- c(0.041, 0.038, 0.042)
weight$dummy <- "Hack"
weight_summary$dummy <- "Hack"
p1 <- ggplot(weight, aes(x=Genotype, y=Fresh_weight, fill=Genotype)) +
  geom_boxplot(width=0.3, outlier.color=NA, alpha=0.7) +
  # geom_jitter(position=position_jitter(width=0), size=0.5, alpha=0.3) +
  geom_text(data=weight_summary, aes(x=Genotype, y=y_pos*1.1, label=label),
            inherit.aes=FALSE, size=8/.pt) +
  # geom_blank(data=upper_limits_Lj, aes(y=Fresh_weight), inherit.aes=FALSE) +
  scale_fill_manual(values=colors) +
  facet_grid(Host~dummy, scales = "free_y", switch = "y")+
  main_theme +
  ylab("Shoot fresh weight/plant (g)") +
  theme(
    legend.position="none",
    plot.title = element_text(size=8, face="bold", hjust=0.5),
    strip.text.y = element_text(size=8, face="bold"),
    strip.text.x = element_text(size=8, color = "white"),
    strip.background.x = element_rect(fill = "white"),
    strip.placement = "outside",
    axis.title.x=element_blank(),
    legend.text = element_text(size=8, colour = "black"),
    legend.title = element_text(size=8, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust=1, size=8),
    axis.title.y=element_text(size=8)
  ) +
  scale_x_discrete(labels=c(
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  )) +
  # scale_y_continuous(expand=expansion(mult=c(0,0)), limits=c(0,0.08)) +
  # scale_y_continuous(expand=expansion(mult=c(0,0))) +
  facetted_pos_scales(
    y = list(
      Host == "Lotus" ~ scale_y_continuous(limits = c(0, 0.05), 
                                           expand = c(0, 0)),
      Host == "Hordeum" ~ scale_y_continuous(limits = c(0, 0.49),
                                             expand = c(0, 0))
    )
  )+
  NULL

p1

ggsave("Plot_shootfw_boxplot.pdf", p1, width = 5, height = 6, unit = "cm")
saveRDS(p1, file = "Plot_shootfw_boxplot.rds")
