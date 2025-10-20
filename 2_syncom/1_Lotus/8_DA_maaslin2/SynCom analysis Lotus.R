## Original script by Ib Thorsgaard Jensen (Aalborg University), modified and extended by Anna Martyn (amartyn@mpipz.mpg.de)

# Clean up.
options(warn=-1)
rm(list=ls())

# Load required packages.
pkg <- c("data.table", "magrittr", "ggplot2", "vegan",
         "RColorBrewer", "Maaslin2", "ComplexHeatmap", "colorRamp2","dplyr","tidyr","ggtext","patchwork","scales")
for(pk in pkg) library(pk, character.only = TRUE)

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load structural zero file.
source("Structural_zeros.R")

# Load metadata, ASV table, and taxonomy file.
design <- read.table("LotusCSSP_LjSC_metadata.txt", header=TRUE, sep="\t")
asv_table <- read.table("feature-table_LotusSYM_LjSC.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, comment.char = "", skip = 1)
taxonomy <- read.table("LjSC_taxonomy.txt", sep="\t", header=TRUE, fill=TRUE)

# Subset for genotypes of interest and only keep matched ASVs in asv table.
design_filtered <- design %>%
  filter(Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")) %>%
  mutate(Compartment = recode(Compartment, "rhizo"="Rhizosphere", "endo"="Root", "nod"="Nodules"))

samples_keep <- design_filtered$SampleID
asv_table_filtered <- asv_table[, colnames(asv_table) %in% samples_keep]
asv_table_matched <- asv_table_filtered[grepl("Lj", rownames(asv_table_filtered)), ]

# Calculate library size for each sample (needed for Structural_zeros3)
design_filtered$library_size <- colSums(asv_table_matched[, samples_keep])

# Split data by compartment.
samples_rhizo <- design_filtered$SampleID[design_filtered$Compartment == "Rhizosphere"]
samples_root  <- design_filtered$SampleID[design_filtered$Compartment == "Root"]

asv_table_rhizo <- asv_table_matched[, samples_rhizo, drop=FALSE]
asv_table_root  <- asv_table_matched[, samples_root, drop=FALSE]

meta_rhizo <- design_filtered %>%
  filter(SampleID %in% samples_rhizo) %>%
  column_to_rownames("SampleID")

meta_root <- design_filtered %>%
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

# Save output.
fwrite(Results_rhizo2, file="DA_SynCom_Lotus_rhizo.csv")
fwrite(Results_root2, file="DA_SynCom_Lotus_root.csv")

# Next, we want to make a plot highlighting the ASVs with differential abundance.

#----------------------------------------
# Compute mean RA for WT in Rhizosphere and Root
#----------------------------------------
asv_table_RA <- sweep(asv_table_matched, 2, colSums(asv_table_matched), "/")
asv_RA_long <- as.data.frame(asv_table_RA) %>%
  rownames_to_column("ASVid") %>%
  pivot_longer(cols=-ASVid, names_to="SampleID", values_to="RA") %>%
  left_join(design_filtered %>% select(SampleID, Compartment, Genotype), by="SampleID")

# Filter for WT
asv_RA_WT <- asv_RA_long %>%
  filter(Genotype=="WT") %>%
  group_by(ASVid, Compartment) %>%
  summarise(mean_RA = mean(RA, na.rm=TRUE), .groups="drop")

# Add taxonomy info
asv_RA_WT <- asv_RA_WT %>%
  left_join(taxonomy %>% select(ASVid, order), by="ASVid") %>%
  mutate(order = ifelse(is.na(order), "Unknown", order))

# Order ASVs by taxonomic order
asv_order_levels <- asv_RA_WT %>%
  distinct(ASVid, order) %>%
  arrange(order) %>%
  pull(ASVid)

asv_RA_WT$ASVid <- factor(asv_RA_WT$ASVid, levels = asv_order_levels)

# Set Compartment factor levels
asv_RA_WT$Compartment <- factor(asv_RA_WT$Compartment, levels = c("Rhizosphere", "Root", "Nodules"))

#----------------------------------------
# Taxonomic color bar
#----------------------------------------
colors <- c(
  "Acidobacteriales"   = "#570861", "Actinomycetales" = "#3e0034",  
  "Bacillales" = "#4b0e5e", "Burkholderiales" = "#645394",
  "Caulobacterales" = "#8e3563", "Chitinophagales" = "#b55385",
  "Chloroflexales" = "#CC99BB", "Corynebacteriales" = "#f6cefc",
  "Flavobacteriales" = "#05294a", "Frankiales" = "#114477",
  "Gaiellales" = "#4477AA", "Gemmatimonadales" = "#77AADD",
  "MB-A2-108" = "#117777", "Micrococcales" = "#44AAAA",
  "Micromonosporales" = "#99D6DD", "Nitrospirales" = "#daf0ee",
  "Pedosphaerales" = "#013220", "Propionibacteriales" = "#117744",
  "Pseudomonadales" = "#88CCAA", "Pseudonocardiales" = "#95bb72",
  "Rhizobiales" = "lightyellow", "S085" = "#774411",
  "Solibacterales" = "#DDAA77", "Sphingomonadales" = "#fdbb6b",
  "Streptomycetales" = "#fed5a4", "Subgroup_7" = "#AA4455",
  "TK10" = "#DD7788", "Xanthomonadales" = "#ffc0cb",
  "Unknown" = "darkgrey", "Other" = "lightgrey"
)

tax_bar <- asv_RA_WT %>%
  distinct(ASVid, order)

p_tax <- ggplot(tax_bar, aes(x=ASVid, y=1, fill=order)) +
  geom_tile() +
  scale_fill_manual(values=colors) +
  theme_void() +
  labs(fill = "Bacterial order") +
  theme(legend.position="right",
        legend.text = element_text(color="black", size=20),
        legend.title = element_text(color="black", size=20),
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

p_bubble <- ggplot(da_all, aes(x = ASV_ID, y = Genotype, fill = DA)) +
  geom_point(shape = 21, size = 5, color = "black") +
  scale_fill_manual(
    values = da_colors,
    labels = c("-1" = "Depleted", "0" = "Non-significant", "1" = "Enriched")
  ) +
  facet_grid(Compartment ~ ., scales = "free_x", switch = "y") +
  labs(fill = "Relative abundance mutant vs. WT") +
  labs(y = "Differencial abundance\nin mutants") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y=element_text(color="black", size=20),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey80", color = "grey50"),
    strip.text.y.left = element_text(color = "black", size = 20, hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 20),
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20),
    panel.spacing = unit(0.4, "lines")
  )

p_bubble <- p_bubble +
  theme(
    axis.text.y = ggtext::element_markdown()  # render x-axis text as markdown
  ) +
  scale_y_discrete(labels = function(x) paste0("*", x, "*"))  # wrap each label in italics

#----------------------------------------
# Mean RA bar plots
#----------------------------------------
p_RA <- ggplot(asv_RA_WT, aes(x=ASVid, y=mean_RA)) +
  geom_bar(stat="identity", fill="grey50") +
  facet_wrap(~Compartment, ncol=1, scales="free", switch = "y") +
  labs(y = "Mean relative\nabundance in WT") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.6)) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(color="black", size=20),
        axis.title.y=element_text(color="black", size=20),
        strip.text = element_text(color="black", size=20),
        strip.placement = "outside",
        panel.spacing = unit(0.4, "lines"))


#----------------------------------------
# Combine plots
#----------------------------------------
# final_plot <- p_RA / p_tax / p_bubble + plot_layout(heights=c(2,0.2,2))
# final_plot

# Extract legends
legend_tax <- cowplot::get_legend(
  p_tax + theme(legend.position = "bottom") +
    guides(
      fill = guide_legend(
        ncol = 3, nrow = 3, 
        title.position = "top",
        title.hjust = 0  # centers title above keys
      )
    )
)

legend_bubble <- cowplot::get_legend(
  p_bubble + theme(legend.position = "bottom") +
    guides(
      fill = guide_legend(
        ncol = 2, nrow = 4, 
        title.position = "top",
        title.hjust = 0  # centers title above keys
      )
    )
)

# Remove individual legends from plots
p_tax_clean <- p_tax + theme(legend.position = "none")
p_bubble_clean <- p_bubble + theme(legend.position = "none")

# Combine plots vertically
main_plot <- p_RA / p_tax_clean / p_bubble_clean + plot_layout(heights = c(0.6, 0.05, 0.6))

# Combine the two legends side by side at the bottom
combined_legend <- cowplot::plot_grid(legend_tax, legend_bubble, ncol = 2, rel_widths = c(0.6, 0.4))

# Final figure: main plot + combined legends
final_plot <- cowplot::plot_grid(main_plot, combined_legend, ncol = 1, rel_heights = c(1, 0.1))

final_plot

# Save final plot.
ggsave("LotusSynCom_DA.pdf", plot = final_plot, width = 21, height = 20)
saveRDS(final_plot, file = "LotusSynCom_DA.rds")


