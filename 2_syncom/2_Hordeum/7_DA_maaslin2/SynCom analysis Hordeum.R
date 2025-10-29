## Original script by Ib Thorsgaard Jensen (Aalborg University), modified and extended by Anna Martyn (amartyn@mpipz.mpg.de)

# Clean up.
options(warn=-1)
rm(list=ls())

# Load required packages.
pkg <- c("data.table", "magrittr", "ggplot2", "vegan",
         "RColorBrewer", "Maaslin2", "ComplexHeatmap", "ggh4x",
         "colorRamp2","dplyr","tidyr","ggtext","patchwork","scales")
for(pk in pkg) library(pk, character.only = TRUE)

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load structural zero file.
source("Structural_zeros.R")

# Load metadata, ASV table, and taxonomy file.
design <- read.table("../1_data/without_input/HordeumSC_metadata.txt",
                     header=TRUE, sep="\t")
asv_table <- read.table("../1_data/without_input/HordeumSC_ASVtable.tsv",
                        sep = "\t", header = TRUE, row.names = 1,
                        check.names = FALSE, comment.char = "", skip = 1)
taxonomy <- read.table("../1_data/without_input/CerealSC_taxonomy_May23.txt",
                       sep="\t", header=TRUE, fill=TRUE)

# Subset for genotypes of interest and only keep matched ASVs in asv table.
asv_table_matched <- asv_table[grepl("_", rownames(asv_table)), , drop = FALSE]

# Calculate library size for each sample (needed for Structural_zeros3)
design$library_size <- colSums(asv_table_matched[, design$SampleID])

# Split data by compartment.
samples_rhizo <- design$SampleID[design$Compartment == "Rhizosphere"]
samples_root  <- design$SampleID[design$Compartment == "Root"]

asv_table_rhizo <- asv_table_matched[, samples_rhizo, drop=FALSE]
asv_table_root  <- asv_table_matched[, samples_root, drop=FALSE]

meta_rhizo <- design %>%
  filter(SampleID %in% samples_rhizo) %>%
  column_to_rownames("SampleID")

meta_root <- design %>%
  filter(SampleID %in% samples_root) %>%
  column_to_rownames("SampleID")

# Set genotype factor levels.
meta_rhizo$Genotype <- factor(meta_rhizo$Genotype, levels=c("WT","symrk","ccamk","nsp1","nsp2"))
meta_root$Genotype  <- factor(meta_root$Genotype, levels=c("WT","symrk","ccamk","nsp1","nsp2"))

# Initialise a results table.
Results_rhizo <- data.frame(matrix(NA, nrow=nrow(asv_table_matched), ncol=8))
rownames(Results_rhizo) <- rownames(asv_table_matched)
colnames(Results_rhizo) <- c(paste("Lfc", c("symrk","ccamk","nsp1","nsp2"), sep="_"),
                             paste("DA", c("symrk","ccamk","nsp1","nsp2"), sep="_"))
Results_root <- Results_rhizo

# Structural zero analysis.
S_rhizo <- Structural_zeros3(asv_table_rhizo, meta_rhizo, group="Genotype",
                             ref="WT", min_reads=20, min_present_reps=2)

S_root <- Structural_zeros3(asv_table_root, meta_root, group="Genotype",
                            ref="WT", min_reads=20, min_present_reps=2)

# Differential abundance analysis using Maaslin2.
## Rhizosphere
M_rhizo <- Maaslin2(input_data = asv_table_rhizo,
                    input_metadata = meta_rhizo,
                    output = "Maaslin2_rhizo",
                    plot_heatmap = FALSE,
                    fixed_effects = "Genotype",
                    min_prevalence = 0.1)

res_dt <- data.table(M_rhizo$results)

for(g in c("symrk","ccamk","nsp1","nsp2")){
  res_g <- res_dt[value == g]
  Results_rhizo[res_g$feature, paste0("Lfc_", g)] <- res_g$coef
  Results_rhizo[res_g$feature, paste0("DA_", g)]  <- (res_g$qval < 0.05) * sign(res_g$coef)
}

# Update results for structural zeros
all_DA_rhizo <- Reduce("union", S_rhizo$struc_zero_DA)
S_sign_rhizo <- S_rhizo$struc_zero_table[all_DA_rhizo,1] - S_rhizo$struc_zero_table[all_DA_rhizo,-1]

for(g in c("symrk","ccamk","nsp1","nsp2")){
  idx <- S_rhizo$struc_zero_DA[[g]]
  Results_rhizo[idx, paste0("DA_", g)]  <- S_sign_rhizo[idx, g]
  Results_rhizo[idx, paste0("Lfc_", g)] <- S_sign_rhizo[idx, g] * Inf
}

## Root
M_root <- Maaslin2(input_data = asv_table_root,
                   input_metadata = meta_root,
                   output = "Maaslin2_root",
                   plot_heatmap = FALSE,
                   fixed_effects = "Genotype",
                   min_prevalence = 0.1)

res_dt <- data.table(M_root$results)

for(g in c("symrk","ccamk","nsp1","nsp2")){
  res_g <- res_dt[value == g]
  Results_root[res_g$feature, paste0("Lfc_", g)] <- res_g$coef
  Results_root[res_g$feature, paste0("DA_", g)]  <- (res_g$qval < 0.05) * sign(res_g$coef)
}

# Update results for structural zeros
all_DA_root <- Reduce("union", S_root$struc_zero_DA)
S_sign_root <- S_root$struc_zero_table[all_DA_root,1] - S_root$struc_zero_table[all_DA_root,-1]

for(g in c("symrk","ccamk","nsp1","nsp2")){
  idx <- S_root$struc_zero_DA[[g]]
  Results_root[idx, paste0("DA_", g)]  <- S_sign_root[idx, g]
  Results_root[idx, paste0("Lfc_", g)] <- S_sign_root[idx, g] * Inf
}

# Merge with taxonomy and relative abundance info.
RA_rhizo <- t(t(asv_table_rhizo)/colSums(asv_table_rhizo))
RA_root  <- t(t(asv_table_root)/colSums(asv_table_root))

Results_rhizo2 <- data.table(ASV_ID = rownames(Results_rhizo), Results_rhizo)
Results_rhizo2 <- merge(Results_rhizo2, taxonomy, by.x="ASV_ID", by.y="ASVid")
Results_rhizo2 <- merge(Results_rhizo2, data.table(ASV_ID = rownames(RA_rhizo), RA_rhizo), by="ASV_ID")

Results_root2 <- data.table(ASV_ID = rownames(Results_root), Results_root)
Results_root2 <- merge(Results_root2, taxonomy, by.x="ASV_ID", by.y="ASVid")
Results_root2 <- merge(Results_root2, data.table(ASV_ID = rownames(RA_root), RA_root), by="ASV_ID")

# ASVs to keep in visulisation
rhizo_any <- apply(Results_rhizo2[,6:9], 1, any)
root_any <- apply(Results_root2[,6:9], 1, any)
rhizo_any[is.na(rhizo_any)] <- F
root_any[is.na(root_any)] <- F
isolate_keep <- Results_rhizo2$ASV_ID[rhizo_any|root_any]

# Save output.
fwrite(Results_rhizo2, file="DA_SynCom_Hordeum_rhizo.csv")
fwrite(Results_root2, file="DA_SynCom_Hordeum_root.csv")

# Next, we want to make a plot highlighting the ASVs with differential abundance.

#----------------------------------------
# Compute mean RA for WT in Rhizosphere and Root
#----------------------------------------
asv_table_RA <- sweep(asv_table_matched, 2, colSums(asv_table_matched), "/")
asv_RA_long <- as.data.frame(asv_table_RA) %>%
  rownames_to_column("ASVid") %>%
  pivot_longer(cols=-ASVid, names_to="SampleID", values_to="RA") %>%
  left_join(design %>% select(SampleID, Compartment, Genotype), by="SampleID")

# Filter for WT
asv_RA_WT <- asv_RA_long %>%
  filter(Genotype=="WT") %>%
  group_by(ASVid, Compartment) %>%
  summarise(mean_RA = mean(RA, na.rm=TRUE), .groups="drop")

# Add taxonomy info
asv_RA_WT <- asv_RA_WT %>%
  left_join(taxonomy %>% select(ASVid, Order), by="ASVid") %>%
  mutate(Order = ifelse(is.na(Order), "Unknown", Order))

# Order ASVs by taxonomic order
asv_order_levels <- asv_RA_WT %>%
  distinct(ASVid, Order) %>%
  arrange(Order) %>%
  pull(ASVid)

asv_RA_WT$ASVid <- factor(asv_RA_WT$ASVid, levels = asv_order_levels)

# Set Compartment factor levels
asv_RA_WT$Compartment <- factor(asv_RA_WT$Compartment, levels = c("Rhizosphere", "Root"))

#----------------------------------------
# Taxonomic color bar
#----------------------------------------
colors <- c(
  "Acidobacteriales"   = "#570861",   # deep purple
  "Actinomycetales"    = "#3e0034",  
  "Bacillales"         = "#4b0e5e",  
  "Burkholderiales"    = "#645394",   # purple-blue
  "Caulobacterales"    = "#8e3563",   # magenta
  "Chitinophagales"    = "#b55385",   # rose
  "Chloroflexales"     = "#CC99BB",   # light purple-pink
  "Corynebacteriales"  = "#f6cefc",
  "Enterobacterales"   = "#191551",
  "Flavobacteriales"   = "#05294a",   # navy
  "Frankiales"         = "#114477",   # dark teal-blue
  "Gaiellales"         = "#4477AA",   # medium blue
  "Gemmatimonadales"   = "#77AADD",   # light blue
  "MB-A2-108"          = "#117777",   # teal
  "Micrococcales"      = "#44AAAA",   # turquoise,
  "Micromonosporales"  = "#99D6DD",   
  "Nitrospirales"      = "#daf0ee",   # pale aqua
  "Pedosphaerales"     = "#013220",   # very dark green
  "Propionibacteriales"= "#117744",   # forest green
  "Pseudomonadales"    = "#88CCAA",   # pastel green
  "Pseudonocardiales"  = "#95bb72",   # lime green (stays in the green cluster)
  "Rhizobiales"        = "#fdbb6b",
  "Rhodobacterales"    = "#C3834D",
  "Rhodospirillales"   = "#302018",
  "S085"               = "#774411",   # brown
  "Solibacterales"     = "#DDAA77",   # beige-brown
  "Sphingobacteriales" = "#8A6642",   
  "Sphingomonadales"   = "lightyellow",
  "Streptomycetales"   = "#fed5a4",   
  "Subgroup_7"         = "#AA4455",   # dark red
  "TK10"               = "#DD7788",   # reddish-pink
  "Xanthomonadales"    = "#ffc0cb",   # light pink
  "Unknown"            = "darkgrey",
  "Other"              = "lightgrey"
)

tax_bar <- asv_RA_WT %>%
  distinct(ASVid, Order) %>% 
  filter(ASVid %in% isolate_keep)

p_tax <- ggplot(tax_bar, aes(x=ASVid, y=1, fill=Order)) +
  geom_tile() +
  scale_fill_manual(values=colors) +
  theme_void() +
  labs(fill = "Bacterial order") +
  theme(legend.position="right",
        legend.text = element_text(color="black", size=8),
        legend.title = element_text(color="black", size=8),
        plot.margin = margin(c(0.5, 0, 0.5, 0), unit = "lines")
        )

#----------------------------------------
# Bubble plot preparation
#----------------------------------------
da_rhizo <- Results_rhizo2 %>%
  select(ASV_ID, DA_symrk, DA_ccamk, DA_nsp1, DA_nsp2) %>%
  pivot_longer(cols=-ASV_ID, names_to="Genotype", values_to="DA") %>%
  mutate(Compartment="Rhizosphere",
         Genotype = recode(Genotype, DA_symrk="symrk", DA_ccamk="ccamk",
                           DA_nsp1="nsp1", DA_nsp2="nsp2"))

da_root <- Results_root2 %>%
  select(ASV_ID, DA_symrk, DA_ccamk, DA_nsp1, DA_nsp2) %>%
  pivot_longer(cols=-ASV_ID, names_to="Genotype", values_to="DA") %>%
  mutate(Compartment="Root",
         Genotype = recode(Genotype, DA_symrk="symrk", DA_ccamk="ccamk",
                           DA_nsp1="nsp1", DA_nsp2="nsp2"))

da_all <- bind_rows(da_rhizo, da_root)
da_all$ASV_ID <- factor(da_all$ASV_ID, levels = asv_order_levels)
da_all$Genotype <- factor(da_all$Genotype, levels = c("symrk","ccamk","nsp1","nsp2"))

da_all$Genotype <- factor(da_all$Genotype, levels = rev(c("symrk","ccamk","nsp1","nsp2")))

da_all <- da_all %>% 
  replace_na(list(DA = 0))

da_all$DA <- factor(da_all$DA, levels = c(-1, 0, 1))
da_colors <- c("-1" = "darkblue", "0" = "white", "1" = "red")

da_all <- da_all %>% filter(ASV_ID %in% isolate_keep)
p_bubble <- ggplot(da_all, aes(x = ASV_ID, y = Genotype, fill = DA)) +
  geom_point(shape = 21, size = 2, color = "black") +
  scale_fill_manual(
    values = da_colors,
    labels = c("-1" = "Depleted", "0" = "Non-significant", "1" = "Enriched")
  ) +
  facet_grid(Compartment ~ ., scales = "free_x", switch = "y") +
  labs(fill = "Relative abundance mutant vs. WT") +
  labs(y = "Differencial abundance\nin mutants") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                               size = 6, color = "black"),
    axis.title.x = element_blank(),
    # axis.title.y=element_text(color="black", size=8),
    axis.title.y=element_blank(),
    strip.placement = "outside",
    legend.position = "bottom",
    # strip.background = element_rect(fill = "grey80", color = "grey50"),
    # strip.text.y.left = element_text(color = "black", size = 8, hjust = 0.5,
    #                                  face = "bold"),
    # axis.text.y = element_text(color = "black", size = 8),
    axis.text.y = element_blank(),
    legend.text = element_text(color = "black", size = 8),
    legend.title = element_text(color = "black", size = 8),
    # strip.background = element_rect(colour = NA),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(c(0.5, 0, 0.5, 0), unit = "lines"),
    panel.spacing = unit(0.4, "lines")
  )

# p_bubble <- p_bubble +
#   theme(
#     axis.text.y = ggtext::element_markdown()  # render x-axis text as markdown
#   ) +
#   scale_y_discrete(labels = function(x) paste0("*", x, "*"))  # wrap each label in italics

#----------------------------------------
# Mean RA bar plots
#----------------------------------------
asv_RA_WT$Compartment <- as.character(asv_RA_WT$Compartment)
asv_RA_WT$Compartment[asv_RA_WT$Compartment == "Rhizosphere"] <- "Rhizo-\nsphere"
asv_RA_WT_nod <- asv_RA_WT %>% 
  filter(Compartment == "Root") %>% 
  mutate(Compartment = "Nod-\nules") %>% 
  mutate(mean_RA = 0)
asv_RA_WT <- rbind(asv_RA_WT, asv_RA_WT_nod)
asv_RA_WT$Compartment <- factor(asv_RA_WT$Compartment, 
                                levels = c("Rhizo-\nsphere", "Root", "Nod-\nules"))
asv_RA_WT <- asv_RA_WT %>% filter(ASVid %in% isolate_keep)
# p_RA <- ggplot(asv_RA_WT, aes(x=ASVid, y=mean_RA)) +
#   geom_bar(stat="identity", fill="grey50") +
#   facet_wrap(~Compartment, ncol=1, strip.position = "left") +
#   labs(y = "Mean relative\nabundance in WT") +
#   # scale_y_continuous(expand = c(0, 0), limits = c(0,0.4)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 0.3)) +
#   theme_bw() +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.y=element_text(color="black", size=8),
#         # axis.title.y=element_text(color="black", size=8),
#         axis.title.y=element_blank(),
#         # strip.text = element_text(color="black", size=8, face = "bold"),
#         strip.placement = "outside",
#         strip.text = element_blank(),
#         strip.background = element_blank(),
#         panel.spacing = unit(0.4, "lines"))+
#   force_panelsizes(cols = c(1, 1, 1), rows = c(1, 1, 0.7))+
#   # facetted_pos_scales(
#   #   y = list(
#   #     Compartment == "Nod-\nules" ~ scale_y_continuous(breaks = c(0, 0.1, 0.2))
#   #   )
#   # )+
#   NULL

# asv_RA_WT <- asv_RA_WT %>% filter(Compartment != "Nod-\nules")
p_RA <- ggplot(asv_RA_WT, aes(x = ASVid, y = mean_RA)) +
  geom_bar(stat="identity", fill="grey50") +
  facet_wrap(~Compartment, ncol=1, scales="free_y",
             strip.position = "left", space = "free_y") +
  labs(y = "Mean relative\nabundance in WT") +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ggtitle("Hordeum")+
  theme(axis.text.x=element_blank(),
        plot.title = element_text(color="black", size=8, face = "bold"),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        # axis.text.y=element_text(color="black", size=8),
        axis.text.y=element_blank(),
        # axis.title.y=element_text(color="black", size=8),
        axis.title.y=element_blank(),
        # strip.text = element_text(color="black", size=8, face = "bold"),
        strip.placement = "outside",
        # strip.background = element_rect(colour = NA),
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = margin(c(0.5, 0, 0.5, 0), unit = "lines"),
        panel.spacing = unit(0.4, "lines"))+
  force_panelsizes(cols = c(1, 1, 1), rows = c(1, 1, 0.7))+
  facetted_pos_scales(
    y = list(
      Compartment == "Nod-\nules" ~ scale_y_continuous(breaks = c(),
                                                       limits = c(0, 0.6),
                                                       expand = c(0, 0)),
      Compartment == "Rhizo-\nsphere" ~ scale_y_continuous(limits = c(0, 0.33),
                                                           expand = c(0, 0)),
      Compartment == "Root" ~ scale_y_continuous(limits = c(0, 0.33),
                                                 expand = c(0, 0))
    )
  )+
  NULL

#----------------------------------------
# Combine plots
#----------------------------------------
# Remove individual legends from plots
p_tax_clean <- p_tax + theme(legend.position = "none")

# Combine plots vertically
main_plot <- p_RA / p_tax_clean / p_bubble + 
  plot_layout(heights = c(0.54, 0.04, 0.42))

main_plot <- main_plot +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing = unit(0, "cm"))

final_plot <- main_plot
final_plot

# Save final plot.
ggsave("HordeumSynCom_DA.pdf", plot = final_plot, 
       width = 21, height = 20, units = "cm")
saveRDS(final_plot, file = "../8_final_figures/HordeumSynCom_DA.rds")

saveRDS(p_RA, file = "../8_final_figures/p_RA_Hv.rds")
saveRDS(p_tax_clean, file = "../8_final_figures/p_tax_clean_Hv.rds")
saveRDS(p_bubble, file = "../8_final_figures/p_bubble_Hv.rds")

