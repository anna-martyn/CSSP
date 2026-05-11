# Seup ------------------------------------------------------------------------
# Cleaning up
options(warn = -1)
rm(list = ls())

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading packages
pkg <- c("ggplot2", "dplyr", "tidyr", "tibble", "ggh4x", "multcompView", "scales")
for(pk in pkg){
  library(pk, character.only = TRUE)
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
  file = "../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv",
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
  file = "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep = "\t",
  header = TRUE,
  fill = TRUE
)

# Cleaning up taxonomy
rename_tax <- function(tax_table){
  colnames(tax_table)[colnames(tax_table) == "Feature.ID"] <- "ASVid"
  tax_table %>%
    separate(
      col = Taxon,
      into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
      sep = "; ",
      fill = "right"
    ) %>%
    mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
    replace(is.na(.), "Unknown") %>%
    select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)
}

lotus_taxonomy <- rename_tax(lotus_taxonomy)
hordeum_taxonomy <- rename_tax(hordeum_taxonomy)

# Confirming that all samples are in both ASV table and metadata
setequal(lotus_design$SampleID, colnames(lotus_asv_table))
setequal(hordeum_design$SampleID, colnames(hordeum_asv_table))

# Ensuring that sampleIDs are in the same order in metadata and ASV table
lotus_asv_table <- lotus_asv_table[,lotus_design$SampleID]
hordeum_asv_table <- hordeum_asv_table[,hordeum_design$SampleID]

# ASV table with relative abundances (RA)
reads_to_RA <- function(asv_table, taxonomy){
  df <- sweep(asv_table, 2, colSums(asv_table), "/") %>%
    # as.data.frame() %>%
    rownames_to_column(var = "ASVid") %>%
    left_join(taxonomy %>% select(ASVid, Order), by = "ASVid")
  df
}

lotus_asv_table_RA <- reads_to_RA(lotus_asv_table, lotus_taxonomy)
hordeum_asv_table_RA <- reads_to_RA(hordeum_asv_table, hordeum_taxonomy)

# Identifying top orders ------------------------------------------------------
# Function to convert ASV table to long form, merge with design, 
# and keep only WT samples
get_long_form_wt <- function(asv_table, design){
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
    filter(Genotype == "WT")

  return(long_form_asv_table)
}

# Function to identify top most abundant orders above a threshold in WT across 
# compartments
get_top_abn_orders <- function(asv_table, design, top, thresh){
  
  # Long form and filter WT
  asv_table_long_form <- get_long_form_wt(asv_table, design)

  # Aggregating ASV table at order-level
  asv_table_order_summary <- asv_table_long_form %>%
    group_by(Order, sampleID, Soil, Compartment) %>%
    summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

  # Identifying top orders with RA above threshold
  top_orders <- asv_table_order_summary %>%
    filter(Compartment %in% c("Rhizosphere", "Root")) %>%
    group_by(Order) %>%
    summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(MeanRA)) %>%
    slice_head(n = top) %>%
    filter(MeanRA >= thresh) %>%
    pull(Order)

  return(top_orders)
}

# Using function to extract top 20 orders above 0.01 RA
top_orders_lotus <- get_top_abn_orders(
  asv_table = lotus_asv_table_RA,
  design = lotus_design,
  top = 20,
  thresh = 0.01
)
top_orders_hordeum <- get_top_abn_orders(
  asv_table = hordeum_asv_table_RA,
  design = hordeum_design,
  top = 20,
  thresh = 0.01
)

# Get orders that are top 20 orders above 0.01 RA in at least one plant
top_orders <- union(top_orders_lotus, top_orders_hordeum)

# Get long form ASV tables
lotus_asv_table_RA_long <- get_long_form_wt(lotus_asv_table_RA, lotus_design)
hordeum_asv_table_RA_long <- get_long_form_wt(hordeum_asv_table_RA, hordeum_design)

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
  group_by(sampleID, Plant, Compartment, Soil, Order) %>%
  summarise(RA = sum(RA), .groups = "drop")

order_table_means <- order_table_RA_long %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  summarise(RA = mean(RA), .groups = "drop") %>%
  mutate(
    Order = factor(Order, levels = c(sort(unique(Order[Order!="Other"])), "Other")),
    Plant = factor(Plant, levels = c("Lotus", "Hordeum")),
    Compartment = factor(Compartment, levels = c("Rhizosphere", "Root", "Nodules")),
    Soil = factor(Soil, levels = c("NPK", "PK", "UF"))
  )


# Main theme for stacked barplot
main_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth=1),
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

# Loading orders colours
colors <- read.table(
  file = "../../../0_files/Bacterial_order_colors.csv",
  header = TRUE,
  sep = ",",
  comment.char = ""
)

# Stacked barplot
bar_plot <- ggplot(order_table_means, aes(x = Soil, y = RA, fill = Order)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand = c(0, 0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill = "Bacterial order") +
  xlab("") +
  guides(fill = guide_legend(nrow = 21)) +
  facet_nested(~ Plant + Compartment, scales = "free_x", space = "free_x") +
  theme(
    axis.text.x = element_text(size = 6, color = "black", angle = 0, vjust = 1),
    strip.text.x = element_text(size = 6, face = "bold"),
    legend.key.size = unit(0.25, 'cm'),
    legend.margin = margin(l = -8)
  )

# Saving plot
ggsave(
  filename = "2_figures/LotusHordeum_Askov_WT_stackedbp_top20_meanRA.pdf",
  plot = bar_plot,
  width = 12,
  height = 6,
  unit = "cm"
)

saveRDS(
  object = bar_plot,
  file = "1_rds_files/LotusHordeum_Askov_WT_stackedbp_top20_meanRA.rds"
)

# Heatmap - soil effects ------------------------------------------------------
# Excluding nodules
order_table_means <- order_table_means %>% filter(Compartment != "Nodules")

# Settting factor levels
non_orders <- c("Unknown", "Other")
order_levels <- unique(order_table_RA_long$Order)
order_levels <- c(setdiff(order_levels, non_orders), non_orders)

order_table_means <- order_table_means %>%
  mutate(
    Order = factor(Order, levels = order_levels),
    Compartment = droplevels(Compartment)
  )

## Hypothesis test ------------------------------------------------------------
letters_list <- list()
opt <- expand.grid(
  Plant = c("Lotus", "Hordeum"),
  Compartment = c("Rhizosphere", "Root"),
  Order = order_levels
)
j <- 1
# Looping over each plant-compartment-order and performing ANOVA
for(i in 1:nrow(opt)){
  # Setting current plant, compartment, and order
  current_plant <- opt$Plant[i]
  current_compartment <- opt$Compartment[i]
  current_order <- opt$Order[i]

  # Keeping only the relavant data
  df_samples <- order_table_RA_long %>% 
    filter(
      Plant == current_plant, 
      Compartment == current_compartment,
      Order == current_order
    )
  
  # Setting UF as the reference soil
  df_samples$Soil <- factor(df_samples$Soil, levels = c("UF", "NPK", "PK"))
  
  # ANOVA
  a <- aov(RA ~ Soil, data = df_samples)
  anv <- anova(a)
  p_val <- anv$`Pr(>F)`[1]

  # Perofrming pair-wise tests with Tukey HSD If global soil effect 
  # is statistically significant, and identifying letters
  if(p_val < 0.05){
    tk <- TukeyHSD(a)
    mcletters <- multcompLetters(tk$Soil[,"p adj"])
    letters <- mcletters$Letters
    df_res <- data.frame(
      Plant = current_plant,
      Compartment = current_compartment,
      Order = current_order,
      Soil = names(letters),
      letter = letters,
      stringsAsFactors = FALSE
    )
    letters_list[[j]] <- df_res
    j <- j+1
  }
}

# Combining all letters into a dataframe
df_letters <- bind_rows(letters_list)

# Joining  significance letters order-level data
order_table_means <- order_table_means %>%
  left_join(df_letters, by = c("Plant", "Compartment", "Order", "Soil"))

## Constructing figure --------------------------------------------------------
# Defining axis breaks and colours for heatmap
breaks <- c(0, 0.005, 0.052, 0.052001, 0.15999, 0.16, 0.34, 0.64)
colors <- c(
  "#1F78B4", "#A6CEE3", "white","#FFFF99",
  "#FF7F00", "#FB9A99", "#E31A1C", "#902121"
)

# Rescale breaks to 0-1 for gradient.
values <- rescale(breaks, to = c(0,1))

# Identifying orders where no pair-wise soil effects were detected
letter_keep <- order_table_means %>% 
  group_by(Plant, Compartment, Order) %>% 
  summarise(remove = all(letter == "a")) %>%
  ungroup() %>%
  mutate(remove = ifelse(is.na(remove), F, remove))

# Removing letters from orders where no pair-wise soil effect were detected
order_table_means <- order_table_means %>% 
  left_join(letter_keep) %>%
  mutate(letter = ifelse(remove, NA, letter))

# Removing "Uknown" and "Other" - these will not be included in heatmap
order_table_means <- order_table_means %>% 
  filter(!(Order %in% c("Other", "Unknown"))) %>% 
  mutate(Order = droplevels(Order))

order_table_means <- order_table_means %>%
  mutate(
    Compartment = ifelse(Compartment == "Rhizosphere", "Rhizo-\nsphere", "Root")
  ) %>%
  mutate(
    Compartment = factor(Compartment, levels = c("Rhizo-\nsphere", "Root"))
  )

# Heatmap
heat_map <- ggplot(order_table_means, aes(x = Soil, y = Order, fill = RA)) +
  geom_tile(color = "grey50") +
  geom_text(aes(label = letter), na.rm = TRUE, size = 6 / .pt) +
  scale_fill_gradientn(
    colors = colors,
    values = values,
    limits = c(0, max(order_table_means$RA, na.rm = TRUE)),
    name = "Relative abundance"
  ) +
  scale_y_discrete(
    limits = rev(levels(order_table_means$Order)),
    position = "right"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "right",
      barwidth = 6,
      barheight = 0.5
    )
  ) +
  main_theme +
  facet_nested(
    ~ Plant + Compartment,
    scales = "free_x",
    space = "free",
    strip = strip_nested(size = "variable")
  ) +
  xlab(NULL) +
  ylab("Bacterial order") +
  theme(
    axis.text.x = element_text(
      size = 6,
      angle = 90,
      vjust = 1,
      hjust = 0.5,
      colour = "black"
    ),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 6, colour = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 6, colour = "black"),
    legend.position = "bottom",
    plot.margin = margin(r = 2.5, l = 15),
    panel.spacing.x = unit(0.1, "lines")
  )

# Saving figure
ggsave(
  filename = "2_figures/LotusHordeum_Askov_WT_orders_heatmap.pdf",
  plot = heat_map,
  width = 12,
  height = 9,
  unit = "cm"
)
saveRDS(
  object = heat_map,
  file = "1_rds_files/LotusHordeum_Askov_WT_orders_heatmap.rds"
)
