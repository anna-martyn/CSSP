# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data.# Clean up.
options(warn=-1)
design <- read.table("BarleyCSSP_SConly_metadata_NEW.txt", header=TRUE, sep="\t")
asv_table <- read.table("feature-table_BarleyCSSP_CerealSConly.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, comment.char = "", skip = 1)
taxonomy <- read.table("CerealSC_taxonomy_new_May23.txt", sep="\t", header=TRUE, fill=TRUE)

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

# Modify design and asv table file to only keep genotypes of interest, and to only keep matched ASVs.
design_filtered <- design %>%
  filter(Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")) %>%
  mutate(Compartment = recode(Compartment, "rhizo"="Rhizosphere", "endo"="Root"))
samples_keep <- design_filtered$SampleID
asv_table_filtered <- asv_table[, colnames(asv_table) %in% samples_keep]

asv_table_matched <- asv_table_filtered[grepl("_", rownames(asv_table_filtered)), , drop = FALSE]

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
  left_join(taxonomy %>% select(ASVid, Order), by="ASVid")

# Reshape to long format and add design info.
df_long <- df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="SampleID", values_to="RA") %>%
  left_join(design_filtered %>% select(SampleID, Compartment, Genotype), by="SampleID")

# Summarise mean relative abundance for wach order per genotype-compartment combination.
df_summary <- df_long %>%
  group_by(Order, Genotype, Compartment) %>%
  summarise(mean_RA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
  mutate(Order = ifelse(is.na(Order), "Unclassified", Order)) %>%
  group_by(Genotype, Compartment) %>%
  mutate(mean_RA = mean_RA / sum(mean_RA)) %>%
  ungroup()

# Remove orders that have relative abundance of 0 in all samples.
df_summary <- df_summary %>%
  group_by(Order) %>%
  filter(sum(mean_RA) > 0) %>%
  ungroup()

# Set factor level for compartments.
df_summary <- df_summary %>%
  mutate(Compartment = factor(Compartment, levels=c("Rhizosphere","Root")))

# Define colours for orders.
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

# Set genotype order and make mutant names italic.
df_summary$Genotype <- factor(df_summary$Genotype,
                              levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2"))

genotype_labels <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

# Plot stacked barplot and save.
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

p1 <- ggplot(df_summary, aes(x=Genotype, y=mean_RA, fill=Order)) +
  geom_bar(stat="identity", width=0.5) +
  facet_wrap(~Compartment, scales="free_x") +
  scale_fill_manual(values=colors) +
  
  scale_x_discrete(labels = genotype_labels) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") + 
  theme(axis.text.x = element_markdown(size=8, color="black", angle=30, hjust=1),
        strip.text.x=element_text(size=8, face="bold"),
        axis.title.x=element_blank(),
        legend.position = "bottom",
        legend.title.position = "top",
        legend.key.size = unit(0.25, "cm"),
        legend.key.spacing.y = unit(0, 'cm'),
        legend.text = element_text(margin = margin(l = 1)),
        # legend.justification = "left",
        legend.justification = c(0.75, 0),     # anchor at bottom left
        legend.direction = "horizontal") +
  guides(fill=guide_legend(ncol = 3))+
  NULL

p1

ggsave("HordeumSC_order_RA_stackedbp.pdf", p1, width=12, height=6, units = "cm")
saveRDS(p1, file="HordeumSC_order_RA_stackedbp.rds")
saveRDS(p1, file="../8_final_figures/HordeumSC_order_RA_stackedbp.rds")

# Next we want to plot the relative abundances of all orders in the different compartment-genotype combinations using barplots.

# Again, look at the relative abundance of each order in each sample.
df_order_sample <- df_long %>%
  mutate(Order = ifelse(is.na(Order), "Unclassified", Order)) %>%
  group_by(Order, SampleID, Compartment, Genotype) %>%
  summarise(RA = sum(RA), .groups = "drop")

# Remove orders that are zero in all samples.
orders_nonzero <- df_order_sample %>%
  group_by(Order) %>%
  summarise(total_RA = sum(RA), .groups = "drop") %>%
  filter(total_RA > 0) %>%
  pull(Order)

df_order_sample <- df_order_sample %>%
  filter(Order %in% orders_nonzero)

# Set the genotype factor levels.
df_order_sample$Genotype <- factor(df_order_sample$Genotype,
                                   levels = c("WT","symrk","ccamk","nsp1","nsp2"))

# Perform significance analysis (ANOVA and Tukey HSD) to look at differences in RA of orders among genotypes in each compartment.

final_results_list <- list()
tukey_results_list <- list()

compartments <- unique(df_order_sample$Compartment)
orders <- unique(df_order_sample$Order)

for(comp in compartments){
  for(ord in orders){
    df_sub <- df_order_sample %>%
      filter(Compartment == comp, Order == ord)

    # Skip if less than 2 genotypes with data
    genotypes_with_data <- df_sub %>%
      group_by(Genotype) %>%
      summarise(n_samples = n(), .groups="drop") %>%
      filter(n_samples > 0) %>%
      pull(Genotype)

    if(length(genotypes_with_data) < 2){
      message("Skipping order ", ord, " in compartment ", comp,
              " because less than 2 genotypes with data")
      next
    }

    # ANOVA.
    ano <- aov(RA ~ Genotype, data=df_sub)
    p_value <- summary(ano)[[1]][["Pr(>F)"]][1]

    # Tukey HSD.
    tukey <- TukeyHSD(ano)

    # Check if any p-values are NA.
    if(any(is.na(tukey$Genotype[,4]))){
      message("Skipping multcompLetters for order ", ord, " in compartment ", comp,
              " because Tukey HSD contains NA p-values")
      tukey_letters <- setNames(rep(NA, 5), c("WT","symrk","ccamk","nsp1","nsp2"))
    } else {
      tukey_letters <- multcompLetters(tukey$Genotype[,4])$Letters
    }

    # Save Tukey full results for this order/compartment.
    tukey_df <- as.data.frame(tukey$Genotype) %>%
      rownames_to_column("Comparison") %>%
      mutate(Order=ord, Compartment=comp)
    tukey_results_list[[paste(ord,comp,sep="_")]] <- tukey_df

    # Save letters for plotting.
    result_row <- data.frame(
      Order = ord,
      Compartment = comp,
      P_Value = p_value,
      WT     = tukey_letters["WT"]     %||% NA,
      symrk  = tukey_letters["symrk"]  %||% NA,
      ccamk  = tukey_letters["ccamk"]  %||% NA,
      nsp1   = tukey_letters["nsp1"]   %||% NA,
      nsp2   = tukey_letters["nsp2"]   %||% NA,
      stringsAsFactors = FALSE
    )
    final_results_list[[paste(ord,comp,sep="_")]] <- result_row
  }
}

final_results <- bind_rows(final_results_list)
tukey_results <- bind_rows(tukey_results_list)

# Save significance analysis results.
write.csv(final_results, "Hordeum_order_RA_ANOVATukey_letters.csv", row.names=FALSE)
write.csv(tukey_results, "Hordeum_order_RA_TukeyHSD_full.csv", row.names=FALSE)

# Prepare data for plotting.

## Summary for plotting mean and SE (per order-compartment-genotype).
df_order_summary <- df_order_sample %>%
  group_by(Order, Compartment, Genotype) %>%
  summarise(Mean_RA = mean(RA), SE_RA = sd(RA)/sqrt(n()), .groups="drop")

## Prepare letters for plotting.
df_plot_letters <- final_results %>%
  pivot_longer(cols=WT:nsp2, names_to="Genotype", values_to="Letter") %>%
  left_join(df_order_summary, by=c("Order"="Order","Compartment","Genotype"))

## Set colours for genotypes.
colors_geno <- c(
  "WT"     = "#A9C289",
  "symrk"  = "#FEDA8B",
  "ccamk"  = "#FDB366",
  "nsp1"   = "#C0E4EF",
  "nsp2"   = "#6EA6CD"
)

## Set main theme.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size=8, color="black"),
  legend.text=element_text(size=8),
  legend.key=element_blank(),
  axis.title.y=element_text(size=8),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size=8, color="black")
)

# Now make plots.

# 1) Plotting all orders (with non-significant and significant differences).

## Make sure factor levels match
df_order_summary$Order <- factor(df_order_summary$Order, levels = unique(df_order_summary$Order))
df_plot_letters$Order <- factor(df_plot_letters$Order, levels = levels(df_order_summary$Order))
df_plot_letters$Genotype <- factor(df_plot_letters$Genotype, levels = levels(df_order_summary$Genotype))

## Define italic labels for legend.
genotype_labels_legend <- c(
  "WT"     = "WT",
  "symrk"  = "*symrk*",
  "ccamk"  = "*ccamk*",
  "nsp1"   = "*nsp1*",
  "nsp2"   = "*nsp2*"
)

## Compute max RA per order × compartment for y-position.
asterisk_df <- df_order_summary %>%
  group_by(Order, Compartment) %>%
  summarise(
    y_position = max(Mean_RA + SE_RA, na.rm=TRUE) + 0.02,
    .groups="drop"
  ) %>%
  left_join(
    final_results %>% select(Order, Compartment, P_Value),
    by = c("Order"="Order", "Compartment")
  ) %>%
  mutate(
    asterisk = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(asterisk))

## Define dodge object and plot.
dodge <- position_dodge(width = 0.9)

p_all <- ggplot(df_order_summary, aes(x=Order, y=Mean_RA, fill=Genotype)) +
  geom_bar(stat="identity", position=dodge, width=0.8, alpha=0.9) +
  geom_errorbar(aes(ymin=Mean_RA-SE_RA, ymax=Mean_RA+SE_RA),
                width=0.3, position=dodge) +
  geom_text(
    data=df_plot_letters,
    aes(x=Order, y=Mean_RA + SE_RA + 0.01, label=Letter, fill=Genotype),
    position=dodge,
    inherit.aes=FALSE,
    size=5
  ) +
  geom_text(
    data = asterisk_df,
    aes(x = Order, y = y_position, label = asterisk),
    inherit.aes = FALSE,
    size = 6
  ) +
  facet_wrap(~Compartment, scales="free_x") +
  scale_fill_manual(values=colors_geno, labels = genotype_labels_legend) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.7)) +
  labs(x="", y="Relative Abundance", fill="Genotype") +
  main_theme +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.text = element_markdown(),   # italic genotypes
    strip.text = element_text(face="bold", size=rel(1)),  # italic facet labels
    plot.title = element_text(size=8)
  )

p_all

ggsave("Hordeum_order_RA_all_orders.pdf", p_all, width=14, height=6)
saveRDS(p_all, file="Hordeum_order_RA_all_orders.rds")

# Now only plot significant orders.

## Get significant orders and sort alphabetically.
sig_orders <- sort(final_results %>% filter(P_Value < 0.05) %>% pull(Order) %>% unique())

## Filter summary and letters.
df_order_summary_sig <- df_order_summary %>% filter(Order %in% sig_orders)
df_plot_letters_sig <- df_plot_letters %>% filter(Order %in% sig_orders)

## Reset factor levels to ensure alphabetical order.
df_order_summary_sig$Order <- factor(df_order_summary_sig$Order, levels = sig_orders)
df_plot_letters_sig$Order <- factor(df_plot_letters_sig$Order, levels = sig_orders)

## Recalculate y-positions for letters.
df_plot_letters_sig <- df_plot_letters_sig %>%
  mutate(y_pos = Mean_RA + SE_RA + 0.01)

## Recalculate asterisk positions.
asterisk_df_sig <- df_order_summary_sig %>%
  group_by(Order, Compartment) %>%
  summarise(y_position = max(Mean_RA + SE_RA, na.rm = TRUE) + 0.02, .groups="drop") %>%
  left_join(final_results %>% select(Order, Compartment, P_Value), 
            by = c("Order" = "Order", "Compartment")) %>%
  mutate(
    asterisk = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(asterisk))

# Removing non-significant orders
df_order_summary_sig <- df_order_summary_sig %>%
  filter(!(Compartment == "Rhizosphere" & Order == "Sphingomonadales")) %>%
  filter(!(Compartment == "Root" & Order == "Enterobacterales")) %>%
  filter(!(Compartment == "Root" & Order == "Pseudomonadales")) %>%
  mutate(Order = as.character(Order))

df_plot_letters_sig <- df_plot_letters_sig %>%
  filter(!(Compartment == "Rhizosphere" & Order == "Sphingomonadales")) %>%
  filter(!(Compartment == "Root" & Order == "Enterobacterales")) %>%
  filter(!(Compartment == "Root" & Order == "Pseudomonadales")) %>%
  mutate(Order = as.character(Order))

## Plot.
p_sig <- ggplot(df_order_summary_sig, aes(x=Order, y=Mean_RA, fill=Genotype)) +
  geom_bar(stat="identity", position=dodge, width=0.8, alpha=0.9) +
  geom_errorbar(aes(ymin=Mean_RA-SE_RA, ymax=Mean_RA+SE_RA),
                width=0.3, position=dodge) +
  geom_text(
    data=df_plot_letters_sig,
    aes(x=Order, y=y_pos+0.05, label=Letter, fill=Genotype),
    position=dodge,
    inherit.aes=FALSE,
    size=2
  ) +
  # geom_text(
  #   data=asterisk_df_sig,
  #   aes(x=Order, y=y_position, label=asterisk),
  #   inherit.aes=FALSE,
  #   size=6
  # ) +
  facet_wrap(~Compartment, scales="free_x") +
  scale_fill_manual(values=colors_geno, labels=genotype_labels_legend) +
  scale_y_continuous(expand=c(0,0), limits=c(0,0.7)) +
  labs(x="", y="Relative Abundance", fill="Genotype") +
  main_theme +
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.text = element_markdown(),
    plot.title = element_text(size=8),
    strip.text = element_text(size=8, face="bold"),
    legend.position = "none"
  )

p_sig

ggsave("Hordeum_order_RA_sign_orders.pdf", p_sig, 
       width=12, height=6, units = "cm")
saveRDS(p_sig, file="Hordeum_order_RA_sign_orders.rds")
saveRDS(p_sig, file="../8_final_figures/Hordeum_order_RA_sign_orders.rds")

# Comment: Order Enterobacterales has letter 'a' for all genotypes, but shows asterisk '*'.
# This is due to the fact that at pairwise comparison levels, the differences are not significant,
# but at an overall ANOVA levels that looks at the whole variance p<0.05.
