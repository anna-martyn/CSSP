# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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

# Load data.
design <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", header=T, sep="\t")
taxonomy <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv", sep="\t", header=T, fill=T)
asv_table <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  skip = 1,
  comment.char = ""
)

# Clean-up the taxonomy file layout.
taxonomy <- taxonomy %>% rename(ASVid = Feature.ID)
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ", fill = "right") %>%
  mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
  replace(is.na(.), "Unknown") %>%
  select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

# Convert the ASV reads to relative abundances and save as new dataframe.
asv_table_norm <- sweep(asv_table, 2, colSums(asv_table), "/")

# Subset the data by soil types NPK, PK, and UF.
# Then check which ASVs in each dataset are present in at least 2 out of 3 biological replicates (RA>0).
get_soil_asv <- function(soil_type) {
  samples <- design %>% filter(Genotype == "Soil", Soil == soil_type) %>% pull(SampleID)
  asv <- asv_table[, samples]
  asv[rowSums(asv > 0) >= 2, ]
}

asv_NPK <- get_soil_asv("NPK")
asv_PK  <- get_soil_asv("PK")
asv_UF  <- get_soil_asv("UF")

# Make a Venn diagram showing the overlap of present ASVs among soils.
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
  theme(legend.position = "none",
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 2, unit = "lines"))+
  NULL

p

# Save the plot.
ggsave(paste("Hordeum_bulk_Venn_ASVs.pdf", sep=""), p, width=5, height=5, bg="white")
saveRDS(p, file = "Hordeum_bulk_Venn_ASVs.rds")
saveRDS(p, file = "../5_final_figure/Hordeum_bulk_Venn_ASVs.rds")

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

# Combine all this information into a dataframe and save it.
assign_category <- function(asv) {
  if(asv %in% ASVs_only_NPK) return("NPK only")
  if(asv %in% ASVs_only_PK)  return("PK only")
  if(asv %in% ASVs_only_UF)  return("UF only")
  if(asv %in% ASVs_NPK_PK)   return("NPK & PK")
  if(asv %in% ASVs_NPK_UF)   return("NPK & UF")
  if(asv %in% ASVs_PK_UF)    return("PK & UF")
  if(asv %in% ASVs_all3)     return("All soils")
  return(NA)
}

all_ASVs <- unique(c(
  ASVs_only_NPK, ASVs_only_PK, ASVs_only_UF,
  ASVs_NPK_PK, ASVs_NPK_UF, ASVs_PK_UF, ASVs_all3
))

ASV_overlap_df <- taxonomy %>%
  filter(ASVid %in% all_ASVs) %>%
  mutate(Category = sapply(ASVid, assign_category))

ASV_overlap_df$Category <- factor(
  ASV_overlap_df$Category,
  levels = c("NPK only", "PK only", "UF only",
             "NPK & PK", "NPK & UF", "PK & UF",
             "All soils")
)

ASV_overlap_df <- ASV_overlap_df %>%
  arrange(Category)

write.csv(ASV_overlap_df, "HordeumCSSP_bulk_ASV_overlap.csv", row.names = FALSE)
write.csv(ASV_overlap_df, "../6_suppl_files/HordeumCSSP_bulk_ASV_overlap.csv", row.names = FALSE)

# Calculate the top20 bacterial orders among all soil samples.
df_order <- data.frame(
  ASVid = rownames(asv_table_norm),
  Order = taxonomy$Order[match(rownames(asv_table_norm), taxonomy$ASVid)],
  asv_table_norm
)

df.long <- df_order %>%
  pivot_longer(cols = -c(ASVid, Order), names_to = "sampleID", values_to = "RA") %>%
  left_join(design %>% select(SampleID, Genotype), by = c("sampleID" = "SampleID")) %>%
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

# Make a function to get the counts for each bacterial order.
get_order_counts <- function(df) {
  counts <- table(df$Order)
  as.data.frame.table(counts) %>% setNames(c("Order", "Count"))
}

# Define the colour for each bacterial order.
colors <- c(
  "Acidobacteriales"="#570861","Burkholderiales"="#645394","Chloroflexales"="#CC99BB",
  "Corynebacteriales"="#f6cefc","Frankiales"="#114477","Gaiellales"="#4477AA",
  "Gemmatimonadales"="#77AADD","MB-A2-108"="#117777","Micrococcales"="#44AAAA",
  "Nitrospirales"="#daf0ee","Pedosphaerales"="#013220","Propionibacteriales"="#117744",
  "Pseudomonadales"="#88CCAA","Rhizobiales"="#fdbb6b","S085"="#774411",
  "Solibacterales"="#DDAA77","Sphingomonadales"="lightyellow","Subgroup_7"="#AA4455",
  "TK10"="#DD7788","Xanthomonadales"="#ffc0cb","Other"="lightgrey"
)

# Write the function for making pie charts.
plot_pie <- function(count_df, title) {
  ggplot(count_df, aes(x = "", y = Count, fill = Order)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors) +
    theme_void() +
    labs(title = title)
}

# Make a pie chart for each dataframe and save as pdf.
ASV_names <- c(
  "ASVs_only_NPK_df", "ASVs_only_PK_df", "ASVs_only_UF_df",
  "ASVs_NPK_UF_df", "ASVs_PK_UF_df", "ASVs_NPK_PK_df", "ASVs_all3_df"
)

for(name in ASV_names) {
  df <- get(name)
  
  # Add top20 bacterial orders or "Other".
  df$Order <- ifelse(df$Order %in% top20_orders, df$Order, "Other")
  
  # Set the factor levels: alphabetically for all except "Other" last.
  df$Order <- factor(df$Order, levels = c(sort(setdiff(unique(df$Order), "Other")), "Other"))
  
  df_counts <- get_order_counts(df)
  
  # Ensure the counts are in the same factor order.
  df_counts$Order <- factor(df_counts$Order, levels = levels(df$Order))
  
  p <- plot_pie(df_counts, title = name)
  print(p)
  ggsave(paste0(name, "_pie.pdf"), p, width = 5, height = 5, bg = "white")
}

# Now combine all pie charts in one figure. 
plot_titles <- c(
  "NPK only", "PK only", "UF only",
  "NPK & PK", "NPK & UF", "PK & UF", "All soils"
)

plot_list <- list()
for(i in seq_along(ASV_names)) {
  df <- get(ASV_names[i])
  
  # Assign bacterial orders, ensure factor levels match the colours.
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

# Add spacers to create a staggered layout (top row centered over bottom row gaps).
s <- plot_spacer()

# Define the top row with 3 pies, and add spacers at the start and end.
top_row <- s + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + s + 
  plot_layout(ncol = 5, widths = c(0.5,1,1,1,0.5))

# Define the bottom row with 4 pies.
bottom_row <- plot_list[[4]] | plot_list[[5]] | plot_list[[6]] | plot_list[[7]]

# Combine both rows with a single legend and update the text sizes.
final_plot <- top_row / bottom_row + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "right",
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 6)
  )

# Update the pie chart title size and reduce the distance to the pies for both rows.
final_plot <- final_plot & 
  theme(
    plot.title = element_text(
      size = 8, 
      hjust = 0.5,
      margin = margin(t=10, b = 1)
    )
  )

# Display the final plot.
final_plot

# Save the final plot.
ggsave(paste("HordeumCSSP_bulk_ASV_overlap_piecharts.pdf", sep=""), final_plot, width=12, height=6, bg="white", units = "cm")
saveRDS(final_plot, file = "HordeumCSSP_bulk_ASV_overlap_piecharts.rds")
saveRDS(final_plot, file = "../5_final_figure/HordeumCSSP_bulk_ASV_overlap_piecharts.rds")

# Create the same final plot without a legend and save it.
final_plot_nolegend <- final_plot & theme(legend.position = "none")
final_plot_nolegend

ggsave("HordeumCSSP_bulk_ASV_overlap_piecharts_nolegend.pdf", final_plot_nolegend, width = 12, height = 6, bg = "white", units = "cm")
saveRDS(final_plot_nolegend, file = "HordeumCSSP_bulk_ASV_overlap_piecharts_nolegend.rds")
saveRDS(final_plot_nolegend, file = "../5_final_figure/HordeumCSSP_bulk_ASV_overlap_piecharts_nolegend.rds")
