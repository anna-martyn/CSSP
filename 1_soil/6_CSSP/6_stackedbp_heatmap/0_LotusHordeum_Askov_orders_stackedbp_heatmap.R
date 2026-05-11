# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("dplyr", "tidyr", "ggplot2", "tibble", "ggh4x", "scales")
for(pk in pkg){
  library(pk, character.only = T)
}

# Loading data
lotus_design <- read.table(
  file = "../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt",
  header = TRUE,
  sep = "\t"
)
lotus_asv_table <- read.table(
  file = "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = ""
)
lotus_taxonomy <- read.table(
  "../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

hordeum_design <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
  header = TRUE,
  sep = "\t"
)
hordeum_asv_table <- read.table(
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  comment.char = "",
  skip = 1
)
hordeum_taxonomy <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

# Cleaning taxonomy
rename_tax <- function(tax_table){
  colnames(tax_table)[colnames(tax_table) == "Feature.ID"] <- "ASVid"
  tax_table %>%
    separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
             sep = "; ", fill = "right") %>%
    mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
    replace(is.na(.), "Unknown") %>%
    select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)
}

lotus_taxonomy <- rename_tax(lotus_taxonomy)
hordeum_taxonomy <- rename_tax(hordeum_taxonomy)

# Relative abundances
reads_to_RA <- function(asv_table, taxonomy){
  df <- sweep(asv_table, 2, colSums(asv_table), "/") %>%
    rownames_to_column(var = "ASVid") %>%
    left_join(taxonomy %>% select(ASVid, Order), by = "ASVid")
  df
}

lotus_asv_table_RA <- reads_to_RA(lotus_asv_table, lotus_taxonomy)
hordeum_asv_table_RA <- reads_to_RA(hordeum_asv_table, hordeum_taxonomy)

# Identifying top orders ------------------------------------------------------
# Function to convert ASV table to long form, merge with design, 
# and keep only WT samples
get_long_form <- function(asv_table, design){
  long_form_asv_table <- asv_table %>%
    pivot_longer(
      cols = -c(ASVid, Order),
      names_to = "sampleID",
      values_to = "RA"
    ) %>%
    left_join(
      design %>% select(SampleID, Plant, Soil, Genotype, Compartment),
      by = c("sampleID" = "SampleID")
    ) %>%
    filter(Compartment %in% c("Rhizosphere", "Root"))

  return(long_form_asv_table)
}

get_top_abn_orders <- function(asv_table, design, top, thresh){
  
  # Long form and filter WT
  asv_table_long_form <- get_long_form(asv_table, design)

  # Aggregating ASV table at order-level
  asv_table_order_summary <- asv_table_long_form %>%
    group_by(Order, sampleID, Soil, Compartment, Genotype) %>%
    summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

  # Identifying top orders with RA above threshold
  top_orders <- asv_table_order_summary %>%
    filter(Genotype == "WT") %>%
    group_by(Order) %>%
    summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MeanRA)) %>%
    slice_head(n = top) %>%
    filter(MeanRA >= thresh) %>%
    pull(Order)

  return(top_orders)
}

# Using function to extract top 20 orders
top_orders_lotus <- get_top_abn_orders(
  asv_table = lotus_asv_table_RA,
  design = lotus_design,
  top = 20,
  thresh = 0
)
top_orders_hordeum <- get_top_abn_orders(
  asv_table = hordeum_asv_table_RA,
  design = hordeum_design,
  top = 20,
  thresh = 0
)

# Get orders belong to the top 20 orders in at least one plant
top_orders <- union(top_orders_lotus, top_orders_hordeum)
saveRDS(top_orders, "../7_DA_analysis/Orders_to_display.rds")

# Get long form ASV tables
lotus_asv_table_RA_long <- get_long_form(lotus_asv_table_RA, lotus_design)
hordeum_asv_table_RA_long <- get_long_form(hordeum_asv_table_RA, hordeum_design)

# Assigning all orders not in top_orders to "Other"
lotus_asv_table_RA_long <- lotus_asv_table_RA_long %>%
  mutate(Order = if_else(Order %in% top_orders, Order, "Other"))

hordeum_asv_table_RA_long <- hordeum_asv_table_RA_long %>%
  mutate(Order = if_else(Order %in% top_orders, Order, "Other"))

# Combining long-form ASV tables from Lotus and Hordeum
asv_table_RA_long <- bind_rows(
  lotus_asv_table_RA_long,
  hordeum_asv_table_RA_long
)

# Stacked barplot -------------------------------------------------------------
# Aggregating long-form ASV table at order-level
order_table_RA_long <- asv_table_RA_long %>%
  group_by(sampleID, Plant, Compartment, Soil, Genotype, Order) %>%
  summarise(RA = sum(RA), .groups = "drop")

# factor levels
non_orders <- c("Unknown", "Other")
order_levels <- sort(setdiff(top_orders, non_orders))
order_levels <- c(order_levels, non_orders)

plant_levels <- c("Lotus", "Hordeum")
compartment_levels <- c("Rhizosphere", "Root", "Nodules")
soil_levels <- c("NPK", "PK", "UF")
genotype_levels <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")

order_table_means <- order_table_RA_long %>%
  group_by(Plant, Compartment, Soil, Genotype, Order) %>%
  summarise(RA = mean(RA), .groups = "drop") %>%
  mutate(
    Order = factor(Order, levels = order_levels),
    Plant = factor(Plant, levels = plant_levels),
    Compartment = factor(Compartment, levels = compartment_levels),
    Soil = factor(Soil, levels = soil_levels),
    Genotype = factor(Genotype, levels = genotype_levels)
  )

# Main theme
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
  axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(size = 6, color = "black"),
  legend.text = element_text(size = 6, color = "black"),
  legend.key = element_blank(),
  axis.title.y = element_text(size = 6),
  legend.position = "right",
  legend.background = element_blank(),
  text = element_text(family = "sans", size = 6, color = "black")
)

# Loading order colours
colors <- read.table(
  "../../../0_files/Bacterial_order_colors.csv",
  header = TRUE,
  sep = ",",
  comment.char = ""
)

# Mutant genotype names in italics
genotype_labels <- c(
  "WT"    = "WT",
  "symrk" = "italic(symrk)",
  "ccamk" = "italic(ccamk)",
  "nsp1"  = "italic(nsp1)",
  "nsp2"  = "italic(nsp2)"
)

# Stacked barplot
bar_plot <- ggplot(order_table_means, aes(x = Genotype, y = RA, fill = Order)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(labels = function(x) parse(text = genotype_labels[x])) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill = "Bacterial order") +
  xlab("") +
  guides(fill = guide_legend(ncol = 1)) +
  facet_nested(
    ~ Plant + Compartment + Soil,
    scales = "free_x",
    space = "free_x"
  ) +
  theme(
    axis.text.x = element_text(
      size = 6,
      color = "black",
      angle = 90,
      vjust = 1,
      hjust = 1
    ),
    strip.text.x = element_text(size = 6, face = "bold"),
    legend.key.size = unit(0.25, 'cm'),
    legend.key.spacing.y = unit(0, 'cm'),
    legend.margin = margin(l = -8)
  )

# Saving plot
ggsave(
  filename = "2_figures/LotusHordeum_Askov_stackedbp_top20_meanRA.pdf",
  plot = bar_plot,
  width = 21,
  height = 9,
  unit = "cm"
)
saveRDS(
  object = bar_plot,
  file = "1_rds_files/LotusHordeum_Askov_stackedbp_top20_meanRA.rds"
)

# Heatmap - soil effects ------------------------------------------------------
# Remvoing 'Other' and 'unknown' orders and setting factor levels
order_table_means <- order_table_RA_long %>%
  filter(!(Order %in% c("Other", "Unknown"))) %>%
  mutate(
    Plant = factor(Plant, levels = plant_levels),
    Compartment = factor(Compartment, levels = compartment_levels),
    Genotype = factor(Genotype, levels = genotype_levels),
    Order = factor(Order, levels = sort(unique(Order)))
  )

# Linear regression to obtain with genotype as covariate, WT as reference level
mutants <- genotype_levels[-1]

opt <- expand.grid(
  Plant = unique(order_table_means$Plant),
  Order = unique(order_table_means$Order),
  Compartment = unique(order_table_means$Compartment),
  Soil = unique(order_table_means$Soil)
)

df_list <- list()
for(i in 1:nrow(opt)){
  current_plant <- opt$Plant[i]
  current_compartment <- opt$Compartment[i]
  current_soil <- opt$Soil[i]
  current_order <- opt$Order[i]

  df <- order_table_means %>%
    filter(
      Plant == current_plant &
        Compartment == current_compartment &
        Soil == current_soil &
        Order == current_order
    )
  p_vals <- coef(summary(lm(RA~Genotype, data = df)))[-1,"Pr(>|t|)"]
  df_list[[i]] <- tibble(
    Plant = opt$Plant[i],
    Compartment = opt$Compartment[i], 
    Soil = opt$Soil[i],
    Order = opt$Order[i],
    Genotype = gsub("Genotype", "", names(p_vals)),
    p.value = p_vals
  )
}
df_pvals <- bind_rows(df_list)

# Correction for multiple testing
## Note: Before this used BH. Changed to bonferroni to be consistent with ANOVA 
# and Tukey used elsewhere. However we really need to change this so it's done properly.
df_pvals <- df_pvals %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  mutate(p.adj = p.adjust(p.value, method = "bonferroni")) %>%
  ungroup() %>%
  mutate(sig = ifelse(p.adj < 0.05, "*", ""))

# Joining results from analysis with order RA table
order_table_means <- order_table_means %>%
  left_join(
    df_pvals %>% select(-p.value),
    by = c("Plant", "Compartment", "Soil", "Order", "Genotype")
  )

# Defining colours and breaks for heatmap
breaks <- c(0, 0.005, 0.052, 0.052001, 0.15999, 0.16, 0.34, 0.64)
heat_colors <- c(
  "#1F78B4",
  "#A6CEE3",
  "white",
  "#FFFF99",
  "#FF7F00",
  "#FB9A99",
  "#E31A1C",
  "#902121"
)
values <- rescale(breaks, to = c(0, 1))

# Mutant genotype names in italics
genotype_labels_heatmap <- c(
  "WT" = "WT",
  "symrk" = "italic(symrk)",
  "ccamk" = "italic(ccamk)",
  "nsp1" = "italic(nsp1)",
  "nsp2" = "italic(nsp2)"
)

# Factor levels
order_table_means$Genotype <- factor(
  order_table_means$Genotype,
  levels = genotype_levels
)

# Heatmap
heat_map <- ggplot(order_table_means, aes(x = Genotype, y = Order, fill = RA)) +
  geom_tile(color = "grey50") +
  geom_text(aes(label = sig), na.rm = TRUE, size = 3) +
  scale_fill_gradientn(
    colors = heat_colors,
    values = values,
    limits = c(0, max(order_table_means$RA, na.rm = TRUE)),
    name = "Relative abundance"
  ) +
  scale_y_discrete(
    limits = rev(levels(order_table_means$Order)),
    position = "right"
  ) +
  scale_x_discrete(labels = function(x) {
    parse(text = genotype_labels_heatmap[x])
  }) +
  facet_nested(
    ~ Plant + Compartment + Soil,
    scales = "free_x",
    space = "free_x"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    axis.text.x = element_text(size = 6, angle = 90, vjust = 1, hjust = 1),
    axis.title.y = element_text(size = 6, color = "black"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 6, colour = "black"),
    legend.position = "bottom",
    plot.margin = margin(r = 10, l = 20)
  ) +
  xlab(NULL) +
  ylab("Bacterial order") +
  guides(fill = guide_colourbar(barheight = 0.5))

# Saving plot
ggsave(
  filename = "2_figures/LotusHordeum_Askov_orders_heatmap.pdf",
  plot = heat_map,
  width = 15,
  height = 9,
  unit = "cm"
)
saveRDS(object = heat_map, file = "1_rds_files/LotusHordeum_Askov_orders_heatmap.rds")
