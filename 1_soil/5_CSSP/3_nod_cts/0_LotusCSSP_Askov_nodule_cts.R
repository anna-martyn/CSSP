# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages. 
pkg <- c("data.table", "ggplot2", "scales", "grid", "car", "multcompView", "ggh4x")

for(pk in pkg){
  library(pk, character.only = T)
}

# Load the nodule count file.
nod <- fread("LotusCSSP_Askov_nodule_cts_all.txt")

# Set the colours for all genotypes.
colors <- data.frame(
  group = c("WT", "symrk","ccamk","nsp1", "nsp2"), 
  color = c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD")
)

# Set the factor levels for soils and genotypes.
nod[,":="(
  Soil_type = factor(Soil_type, levels = c("NPK", "PK", "UF")),
  Genotype = factor(Genotype, levels = colors$group)
)]

# Set the main theme for the plot.
dodge <- position_dodge(width = 0.9)

main_theme <- theme(
  panel.background=element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text.x = element_text(size = 8, angle = 45, vjust = 1, 
                             hjust=1, colour = "black"),
  axis.text.y = element_text(size = 8, colour = "black"),
  legend.background=element_blank(),
  legend.key=element_blank(),
  text=element_text(family="sans")
)

# Make the plot.
p1 <- ggplot(nod, aes(x=Genotype, y=pink, fill=Genotype)) +
  geom_boxplot(width = 0.5, position = dodge, outlier.size = 0.5)+
  scale_fill_manual(values=as.character(colors$color)) +
  facet_wrap(~Soil_type, nrow=1)+
  main_theme +
  ylab("Pink nodule counts/plant")+
  scale_y_continuous(limits = c(0, 15))+
  theme(
    legend.position="none", 
    plot.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 8, colour = "black"),
    strip.text = element_text(size = 8, colour = "black", face = "bold"),
    legend.text = element_text(size = 8, colour = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(size = 8, angle = 45, vjust = 1,
                               hjust=1, colour = "black"),
    legend.key.size = unit(0.5,"cm")
  )+
  scale_x_discrete(labels=c(
    "WT"="WT",
    "symrk"=expression(italic("symrk")),
    "ccamk"=expression(italic("ccamk")),
    "nsp1"=expression(italic("nsp1")),
    "nsp2"=expression(italic("nsp2"))
  ))
p1

ggsave(paste0("LotusCSSP_AskovSoils_nod_cts.pdf"), p1, width=6, height=5, units = "cm")
saveRDS(p1, "LotusCSSP_AskovSoils_nod_cts.rds")
saveRDS(p1, "../8_final_figures/LotusCSSP_AskovSoils_nod_cts.rds")
