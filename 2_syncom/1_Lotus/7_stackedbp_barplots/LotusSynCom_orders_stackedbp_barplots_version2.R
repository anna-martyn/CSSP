# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data.
design <- read.table("LotusCSSP_LjSC_metadata.txt", header=TRUE, sep="\t")
asv_table <- read.table("feature-table_LotusSYM_LjSC.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, comment.char = "", skip = 1)
taxonomy <- read.table("LjSC_taxonomy.txt", sep="\t", header=TRUE, fill=TRUE)

# Load packages.
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggtext)
library(ggforce)
library(ggalluvial)
library(FSA)
library(multcompView)
library(scales)
library(Maaslin2)

# Modify design and asv table file to only keep genotypes of interest, and to only keep matched ASVs.
design_filtered <- design %>%
  filter(Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")) %>%
  mutate(Compartment = recode(Compartment, "rhizo"="Rhizosphere", "endo"="Root", "nod"="Nodules"))
samples_keep <- design_filtered$SampleID
asv_table_filtered <- asv_table[, colnames(asv_table) %in% samples_keep]

asv_table_matched <- asv_table_filtered[grepl("Lj", rownames(asv_table_filtered)), ]

# Check whether all matched ASVs are present in taxonomy file.
missing_asvs <- setdiff(rownames(asv_table_matched), taxonomy$ASVid)
if(length(missing_asvs) == 0){
  message("All matched ASVs are present in taxonomy.")
} else {
  warning(length(missing_asvs), " matched ASVs are missing in taxonomy:")
  print(missing_asvs)
}

# Convert ASV reads to relative abundances and add taxonomic info at order level.
asv_table_norm <- sweep(asv_table_matched, 2, colSums(asv_table_matched), "/")
df <- as.data.frame(asv_table_norm) %>%
  rownames_to_column(var="ASVid") %>%
  left_join(taxonomy %>% select(ASVid, order), by="ASVid")

# Reshape to long format and add design info.
df_long <- df %>%
  pivot_longer(cols=-c(ASVid, order), names_to="SampleID", values_to="RA") %>%
  left_join(design_filtered %>% select(SampleID, Compartment, Genotype), by="SampleID")

# Summarize mean relative abundance per genotype-compartment combination.
# df_summary <- df_long %>%
#   group_by(order, Genotype, Compartment) %>%
#   summarise(mean_RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
#   mutate(order = ifelse(is.na(order), "Unclassified", order)) %>%
#   group_by(Genotype, Compartment) %>%
#   mutate(mean_RA = mean_RA / sum(mean_RA)) %>%
#   ungroup()

df_summary <- df_long %>% 
  group_by(SampleID, order, Genotype, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop") %>% 
  group_by(order, Genotype, Compartment) %>%
  summarise(mean_RA = mean(RA, na.rm = TRUE),
            sd = sd(RA, na.rm = T),
            N = length(order),
            .groups = "drop")
df_summary$se <- df_summary$sd/sqrt(df_summary$N)

# Ensure compartment order.
df_summary <- df_summary %>%
  mutate(Compartment = factor(Compartment, levels=c("Rhizosphere","Root","Nodules")))

# Define colours for orders.
colors <- c(
  "Acidobacteriales"   = "#570861",   # deep purple
  "Actinomycetales"    = "#3e0034",  
  "Bacillales"         = "#4b0e5e",  
  "Burkholderiales"    = "#645394",   # purple-blue
  "Caulobacterales"    = "#8e3563",   # magenta
  "Chitinophagales"    = "#b55385",   # rose
  "Chloroflexales"     = "#CC99BB",   # light purple-pink
  "Corynebacteriales"  = "#f6cefc",   # very light pink
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
  "S085"               = "#774411",   # brown
  "Solibacterales"     = "#DDAA77",   # beige-brown
  "Sphingomonadales"   = "lightyellow",
  "Streptomycetales"    = "#fed5a4",   # pink-magenta (to match other actinobacteria hues)
  "Subgroup_7"         = "#AA4455",   # dark red
  "TK10"               = "#DD7788",   # reddish-pink
  "Xanthomonadales"    = "#ffc0cb",   # light pink
  "Unknown"            = "darkgrey",
  "Other"              = "lightgrey"
)

# Set genotype order and make mutant names italic.
df_summary$Genotype <- factor(df_summary$Genotype,
                              levels = rev(c("WT", "symrk", "ccamk", "nsp1", "nsp2")))

genotype_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Plot stacked barplot.
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

# p1 <- ggplot(df_summary, aes(x=Genotype, y=mean_RA, fill=order)) +
#   geom_bar(stat="identity", width=0.5) +
#   facet_wrap(~Compartment, scales="free_x", space = "free_x") +
#   scale_fill_manual(values=colors) +
#   scale_x_discrete(labels = genotype_labels) +
#   scale_y_continuous(expand=c(0,0)) +
#   main_theme +
#   ylab("Mean relative abundance") +
#   labs(fill="Bacterial order") +
#   theme(axis.text.x = element_markdown(size=8, color="black", angle=30, hjust=1),
#         strip.text.x=element_text(size=8, face="bold"),
#         legend.position = "bottom",
#         legend.title.position = "top",
#         legend.key.size = unit(0.25, "cm"),
#         legend.key.spacing.y = unit(0, 'cm'),
#         axis.title.x=element_blank()) +
#   guides(fill=guide_legend(nrow=3))+
#   NULL

p1 <- ggplot(df_summary, aes(y=Genotype, x=mean_RA, fill=order)) +
  geom_bar(stat="identity", width=0.5) +
  facet_wrap(~Compartment, scales="free_y", space = "free_y", nrow = 3, 
             strip.position = "right") +
  scale_fill_manual(values=colors) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_discrete(labels = genotype_labels) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") +
  theme(axis.text.x = element_markdown(size=8, color="black"),
        # strip.text.x=element_text(size=8, face="bold"),
        strip.text.y = element_text(angle = 0, size=8, face="bold"),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.key.size = unit(0.25, "cm"),
        legend.key.spacing.y = unit(0, 'cm'),
        axis.title.x=element_blank(),
        axis.text.y = element_markdown()) +
  guides(fill=guide_legend(nrow=3))+
  NULL

p1

# Save plot.
ggsave("LotusSC_order_RA_stackedbp.pdf", p1, width=12, height=6, unit = "cm")
saveRDS(p1, file="LotusSC_order_RA_stackedbp.rds")
saveRDS(p1, file="../10_final_figures/LotusSC_order_RA_stackedbp.rds")

# Order level comparisons with Maaslin2 ----
df_summary$Genotype <- factor(df_summary$Genotype,
                              levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2"))


asv_table_df <- as.data.frame(asv_table_matched)
design_df <- data.frame(design[,-1], row.names = design$SampleID)
taxonomy2 <- data.frame(taxonomy[,-1], row.names = taxonomy$ASVid)
Orders <- taxonomy2[rownames(asv_table_df),"order"]

order_table <- rowsum(asv_table_df, Orders)
rhiz_smp <- rownames(design_df)[design_df$Compartment == "rhizo" & 
                                design_df$Genotype != "f6h1"]
root_smp <- rownames(design_df)[design_df$Compartment == "endo" &
                                design_df$Genotype != "f6h1"]
order_table_rhiz <- order_table[,rhiz_smp]
order_table_root <- order_table[,root_smp]
design_rhiz <- design_df[colnames(order_table_rhiz),]
design_root <- design_df[colnames(order_table_root),]
design_rhiz$Genotype <- factor(
  design_rhiz$Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)
design_root$Genotype <- factor(
  design_root$Genotype, levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
)

order_test_rhiz <- Maaslin2(input_data = order_table_rhiz,
                            input_metadata = design_rhiz,
                            fixed_effects = "Genotype",
                            output = "Maslin_output")

order_test_root <- Maaslin2(input_data = order_table_root,
                            input_metadata = design_root,
                            fixed_effects = "Genotype",
                            output = "Maslin_output")

res_rhiz <- order_test_rhiz$results
res_rhiz <- res_rhiz[res_rhiz$qval<0.05,c("feature", "value")]
colnames(res_rhiz) <- c("order", "Genotype")
res_rhiz$Sig <- "*"
res_rhiz$Compartment <- "Rhizosphere"

res_root <- order_test_root$results
res_root <- res_root[res_root$qval<0.05, c("feature", "value")]
colnames(res_root) <- c("order", "Genotype")
res_root$Sig <- "*"
res_root$Compartment <- "Root"

df_ord_rhiz <- df_summary[df_summary$Compartment == "Rhizosphere",]
df_ord_rhiz <- merge(df_ord_rhiz, res_rhiz, all.x = T,
                     by = c("order", "Compartment", "Genotype"))
df_ord_rhiz$Sig[is.na(df_ord_rhiz$Sig)] <- ""
rhiz_sig_ord <- unique(df_ord_rhiz$order[df_ord_rhiz$Sig == "*"])
df_ord_rhiz <- df_ord_rhiz[df_ord_rhiz$order %in% rhiz_sig_ord,]

df_ord_root <- df_summary[df_summary$Compartment == "Root",]
df_ord_root <- merge(df_ord_root, res_root, all.x = T,
                     by = c("order", "Compartment", "Genotype"))
df_ord_root$Sig[is.na(df_ord_root$Sig)] <- ""
root_sig_ord <- unique(df_ord_root$order[df_ord_root$Sig == "*"])
df_ord_root <- df_ord_root[df_ord_root$order %in% root_sig_ord,]

df_ord <- rbind(df_ord_rhiz, df_ord_root)

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

df_ord$Genotype <- factor(df_ord$Genotype, 
                          levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2"))
p_sig <- ggplot(df_ord, aes(x = order, y = mean_RA, fill = Genotype)) +
  geom_bar(stat = "identity", position = dodge, width = 0.8, alpha = 0.9) +
  geom_errorbar(aes(ymin=mean_RA-1.96*se, ymax=mean_RA+1.96*se),
                width = 0.3, position = dodge) +
  geom_text(
    data = df_ord,
    aes(x = order, y=mean_RA+1.96*se+0.005, label=Sig, fill=Genotype),
    position=dodge,
    inherit.aes=FALSE,
    size=4
  ) + 
  facet_wrap(~Compartment, scales="free_x") +
  scale_fill_manual(values = colors_geno, labels=genotype_labels_legend) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x="", y="Relative Abundance", fill="Genotype") +
  main_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust=1),
    legend.text = element_markdown(),
    strip.text = element_text(face="bold", size=rel(1)),
    legend.position = "none",
    plot.title = element_text(size = 8)
  )

p_sig

ggsave("Lotus_order_RA_sign_orders_asterisks.pdf", 
       p_sig, width=12, height=6, unit = "cm")
saveRDS(p_sig, file="Lotus_order_RA_sign_orders_asterisks.rds")
saveRDS(p_sig, 
        file="../10_final_figures/Lotus_order_RA_sign_orders_asterisks.rds")

