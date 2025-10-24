# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data.
design <- read.table("BarleyCSSP_Askov_reseq_metadata.txt", header=T, sep="\t")
asv_table <- read.table("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table("Barley_Askov_Rep_10_4_taxonomy.txt", sep="\t", header=T, fill=T)

# Load required packages.
library(ggVennDiagram)
library(ggtext)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(patchwork)
library(cowplot)

# Clean-up taxonomy file layout.
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^.{3}", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Remove the three Lotus soil samples from the dataset and only keep samples that are in asv table.
design <- design %>%
  filter(!str_detect(Description, "Lj") & Sample_ID %in% colnames(asv_table))

# Subset and reorder ASV table to match filtered design and only keep ASVs present in taxonomy file.
asv_table <- asv_table %>%
  select(all_of(design$Sample_ID)) %>%
  filter(rownames(.) %in% taxonomy$ASVid)

# Convert ASV reads to relative abundances and save as new dataframe.
asv_table_norm <- sweep(asv_table, 2, colSums(asv_table), "/")

# Subset data by soil types NPK, PK, and UF. Then, check which ASVs in each dataset are present in at least 2 out of 3 biological replicates (RA>0).
get_soil_asv <- function(soil_type) {
  samples <- design %>% filter(Genotype == "Soil", Soil == soil_type) %>% pull(Sample_ID)
  asv <- asv_table[, samples]
  asv[rowSums(asv > 0) >= 2, ]
}

asv_NPK <- get_soil_asv("NPK")
asv_PK  <- get_soil_asv("PK")
asv_UF  <- get_soil_asv("UF")

# Make a Venn diagram of the overlap of present ASVs among soils.
venn_data <- list(
  NPK = rownames(asv_NPK),
  PK = rownames(asv_PK),
  UF = rownames(asv_UF)
)

p <- ggVennDiagram(
  venn_data,
  label = "count",
  label_alpha = 0,
  label_size = 8/.pt,
  set_size = 8/.pt,
  edge_size = 0.5
) +
  scale_fill_gradient(low = "white", high = "white") +
  # labs(title = "Unique <i>vs</i>. Shared ASVs") +
  # theme(plot.title = element_markdown(size = 12, hjust = 0.5), legend.position = "none")+
  theme(legend.position = "none",
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 2, unit = "lines"))+
  # geom_segment(aes(x = 9, xend = 11, y = 1, yend = 3))+
  # geom_segment(aes(x = 9, xend = 11, y = -5, yend = -7))+
  NULL

p

ggsave(paste("Venn_ASVs_bulk.pdf", sep=""), p, width=5, height=5, bg="white")
saveRDS(p, file = "Venn_ASVs_bulk.rds")
saveRDS(p, file = "../6_final_figure/Venn_ASVs_bulk.rds")

# Make a new dataframe for each combination of unique or overlapping ASVs.

## Which ASVs are present in NPK, PK, UF?
NPK_ASVs <- rownames(asv_NPK)
PK_ASVs  <- rownames(asv_PK)
UF_ASVs  <- rownames(asv_UF)

## Make new dataframes for unique/overlapping.
ASVs_only_NPK <- setdiff(NPK_ASVs, union(PK_ASVs, UF_ASVs))
ASVs_only_PK  <- setdiff(PK_ASVs,  union(NPK_ASVs, UF_ASVs))
ASVs_only_UF  <- setdiff(UF_ASVs,  union(NPK_ASVs, PK_ASVs))

ASVs_NPK_UF <- intersect(NPK_ASVs, UF_ASVs) %>% setdiff(PK_ASVs)
ASVs_PK_UF  <- intersect(PK_ASVs, UF_ASVs) %>% setdiff(NPK_ASVs)
ASVs_NPK_PK <- intersect(NPK_ASVs, PK_ASVs) %>% setdiff(UF_ASVs)

ASVs_all3 <- Reduce(intersect, list(NPK_ASVs, PK_ASVs, UF_ASVs))

# Add taxonomic info.
ASVs_only_NPK_df <- taxonomy %>% filter(ASVid %in% ASVs_only_NPK)
ASVs_only_PK_df  <- taxonomy %>% filter(ASVid %in% ASVs_only_PK)
ASVs_only_UF_df  <- taxonomy %>% filter(ASVid %in% ASVs_only_UF)
ASVs_NPK_UF_df   <- taxonomy %>% filter(ASVid %in% ASVs_NPK_UF)
ASVs_PK_UF_df    <- taxonomy %>% filter(ASVid %in% ASVs_PK_UF)
ASVs_NPK_PK_df   <- taxonomy %>% filter(ASVid %in% ASVs_NPK_PK)
ASVs_all3_df     <- taxonomy %>% filter(ASVid %in% ASVs_all3)

# Calculate top20 bacterial orders (all soil samples)
df_order <- data.frame(
  ASVid = rownames(asv_table_norm),
  Order = taxonomy$Order[match(rownames(asv_table_norm), taxonomy$ASVid)],
  asv_table_norm
)

df.long <- df_order %>%
  pivot_longer(cols = -c(ASVid, Order), names_to = "sampleID", values_to = "RA") %>%
  left_join(design %>% select(Sample_ID, Genotype), by = c("sampleID" = "Sample_ID")) %>%
  filter(Genotype == "Soil")

df.long_order <- df.long %>%
  group_by(Order, sampleID) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

top20_orders <- df.long_order %>%
  group_by(Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20) %>%
  pull(Order)

# Function to get counts per Order
get_order_counts <- function(df) {
  counts <- table(df$Order)
  as.data.frame.table(counts) %>% setNames(c("Order", "Count"))
}

# Define colours for orders.
colors <- c(
  "Acidobacteriales"="#570861","Burkholderiales"="#645394","Chloroflexales"="#CC99BB",
  "Corynebacteriales"="#f6cefc","Frankiales"="#114477","Gaiellales"="#4477AA",
  "Gemmatimonadales"="#77AADD","MB-A2-108"="#117777","Micrococcales"="#44AAAA",
  "Nitrospirales"="#daf0ee","Pedosphaerales"="#013220","Propionibacteriales"="#117744",
  "Pseudomonadales"="#88CCAA","Rhizobiales"="#fdbb6b","S085"="#774411",
  "Solibacterales"="#DDAA77","Sphingomonadales"="lightyellow","Subgroup_7"="#AA4455",
  "TK10"="#DD7788","Xanthomonadales"="#ffc0cb","Other"="lightgrey"
)

plot_pie <- function(count_df, title) {
  ggplot(count_df, aes(x = "", y = Count, fill = Order)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors) +
    theme_void() +
    labs(title = title)
}

# Process all ASV dataframes for pie charts.
ASV_names <- c(
  "ASVs_only_NPK_df", "ASVs_only_PK_df", "ASVs_only_UF_df",
  "ASVs_NPK_UF_df", "ASVs_PK_UF_df", "ASVs_NPK_PK_df", "ASVs_all3_df"
)

for(name in ASV_names) {
  df <- get(name)
  
  # Add top20 Orders or "Other"
  df$Order <- ifelse(df$Order %in% top20_orders, df$Order, "Other")
  
  # Set factor levels: alphabetically for all except "Other" last
  df$Order <- factor(df$Order, levels = c(sort(setdiff(unique(df$Order), "Other")), "Other"))
  
  df_counts <- get_order_counts(df)
  
  # Ensure the counts are in the same factor order
  df_counts$Order <- factor(df_counts$Order, levels = levels(df$Order))
  
  p <- plot_pie(df_counts, title = name)
  print(p)
  ggsave(paste0(name, "_pie.pdf"), p, width = 5, height = 5, bg = "white")
}

# Now combine all pie charts in one figure. 

# plot_titles <- c(
#   "NPK only", "PK only", "UF only",
#   "NPK & PK", "NPK & UF", "PK & UF", "NPK, PK & UF"
# )

plot_titles <- c(
  "NPK only", "PK only", "UF only",
  "NPK & PK", "NPK & UF", "PK & UF", "All soils"
)

plot_list <- list()
for(i in seq_along(ASV_names)) {
  df <- get(ASV_names[i])
  
  # Assign Orders, ensure factor levels match colors
  df$Order <- ifelse(df$Order %in% names(colors), df$Order, "Other")
  df$Order <- factor(df$Order, levels = names(colors))
  
  df_counts <- get_order_counts(df)
  df_counts$Order <- factor(df_counts$Order, levels = names(colors))
  
  p <- ggplot(df_counts, aes(x = 1, y = Count, fill = Order)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors, guide = guide_legend(ncol = 1)) +
    theme_void() +
    labs(title = plot_titles[i],fill="Bacterial order") +
    theme(plot.title = element_text(hjust = 0.5, size = 6), margin = margin(b = 1))
  
  plot_list[[i]] <- p
}

# Add spacers to create staggered layout (top row centered over bottom row gaps)
s <- plot_spacer()

# Top row with 3 pies, add spacers at start and end
top_row <- s + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + s + 
  plot_layout(ncol = 5, widths = c(0.5,1,1,1,0.5))

# Bottom row with 4 pies
bottom_row <- plot_list[[4]] | plot_list[[5]] | plot_list[[6]] | plot_list[[7]]

# Combine rows with single legend and updated text sizes
final_plot <- top_row / bottom_row + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 6),    # legend text size
    legend.title = element_text(size = 6)     # legend header size
  )

# Update pie chart titles size and reduce distance to pie (both rows)
final_plot <- final_plot & 
  theme(
    plot.title = element_text(
      size = 8, 
      hjust = 0.5,
      margin = margin(t=10, b = 1)  # very small space below title
    )
  )

# Add overall title
final_plot <- final_plot + plot_annotation(
  # title = "Taxonomic information of unique and shared ASVs",
  theme = theme(plot.title = element_text(size = 12, hjust = 0.5))
)

# Display final plot
final_plot

# Save final plot.
ggsave(paste("Barley_bulk_ASV_overlap_piecharts.pdf", sep=""),
       final_plot, width=12, height=6, bg="white", units = "cm")
saveRDS(final_plot, file = "Barley_bulk_ASV_overlap_piecharts.rds")
saveRDS(final_plot, 
        file = "../6_final_figure/Barley_bulk_ASV_overlap_piecharts.rds")

# Create same final plot without a legend.
final_plot_nolegend <- final_plot & theme(legend.position = "none")
final_plot_nolegend

ggsave("Barley_bulk_ASV_overlap_piecharts_nolegend.pdf", final_plot_nolegend,
       width = 12, height = 6, bg = "white", units = "cm")
saveRDS(final_plot_nolegend, 
        file = "Barley_bulk_ASV_overlap_piecharts_nolegend.rds")
saveRDS(final_plot_nolegend, 
        file = "../6_final_figure/Barley_bulk_ASV_overlap_piecharts_nolegend.rds")
