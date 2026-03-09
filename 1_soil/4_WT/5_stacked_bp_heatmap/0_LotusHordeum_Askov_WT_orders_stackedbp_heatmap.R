# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggh4x)
library(FSA)
library(multcompView)
library(scales)

# Load the Lotus and Hordeum input files.
Lotus_design <- read.table("../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", header=T, sep="\t")
Lotus_asv_table <- read.table( "../../1_data/1_Lotus/LotusCSSP_AskovSoils_ASVtable_10_4_nospike.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = "")
Lotus_taxonomy <- read.table("../../1_data/1_Lotus/LotusCSSP_AskovSoils_taxonomy_10_4.tsv", sep="\t", header=TRUE, fill=TRUE)

Hordeum_design <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt", header=T, sep="\t")
Hordeum_asv_table <- read.table(
  "../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_ASVtable_10_4.tsv",
  sep = "\t", header = TRUE, row.names = 1, 
  check.names = FALSE, comment.char = "", skip = 1)
Hordeum_taxonomy <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_taxonomy_10_4.tsv",
  sep="\t", header=TRUE, fill=TRUE)

# Rename the columns "Feature.ID" to "ASVid" in taxonomy files.
# Clean up the taxonomy file layouts.
rename_tax <- function(tax_table){
  colnames(tax_table)[colnames(tax_table) == "Feature.ID"] <- "ASVid"
  tax_table %>%
    separate(Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
             sep = "; ", fill = "right") %>%
    mutate(across(Kingdom:Species, ~sub("^[a-z]__", "", .))) %>%
    replace(is.na(.), "Unknown") %>%
    select(ASVid, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)
}

Lotus_taxonomy <- rename_tax(Lotus_taxonomy)
Hordeum_taxonomy <- rename_tax(Hordeum_taxonomy)

# Re-order the data matrices (only keep samples in design that are found in ASV table, and then filter for ones present in taxonomy file and adjust).
Lotus_design <- Lotus_design %>%
  filter(SampleID %in% colnames(Lotus_asv_table)) 
Lotus_asv_table <- Lotus_asv_table %>% 
  select(all_of(Lotus_design$SampleID)) %>%
  filter(rownames(.) %in% Lotus_taxonomy$ASVid)

Hordeum_design <- Hordeum_design %>%
  filter(SampleID %in% colnames(Hordeum_asv_table)) 
Hordeum_asv_table <- Hordeum_asv_table %>%
  select(all_of(Hordeum_design$SampleID)) %>%
  filter(rownames(.) %in% Hordeum_taxonomy$ASVid)

# Convert ASV reads to relative abundances and save as new dataframe.
asv_to_df <- function(asv_table, taxonomy){
  df <- sweep(asv_table, 2, colSums(asv_table), "/") %>%
    as.data.frame() %>%
    rownames_to_column(var="ASVid") %>%
    left_join(taxonomy %>% select(ASVid, Order), by="ASVid")
  df
}

Lotus_df <- asv_to_df(Lotus_asv_table, Lotus_taxonomy)
Hordeum_df <- asv_to_df(Hordeum_asv_table, Hordeum_taxonomy)

# Reshape to long format and filter for WT samples.
Lotus_df.long <- Lotus_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Lotus_design %>% select(SampleID, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="WT")

Hordeum_df.long <- Hordeum_df %>%
  pivot_longer(cols=-c(ASVid, Order), names_to="sampleID", values_to="RA") %>%
  left_join(Hordeum_design %>% select(SampleID, Plant, Soil, Genotype, Compartment), by=c("sampleID"="SampleID")) %>%
  filter(Genotype=="WT")

# Next we want to choose which bacterial orders to use for representation of the bacterial community structure in a stacked barplot.

## For this, we will first make a summary of the relative abundances by order per sample-soil-compartment group for Lotus and Hordeum.
Lotus_order_summary <- Lotus_df.long %>%
  group_by(Order, sampleID, Soil, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

Hordeum_order_summary <- Hordeum_df.long %>%
  group_by(Order, sampleID, Soil, Compartment) %>%
  summarise(RA = sum(RA, na.rm = TRUE), .groups = "drop")

## Next, we look at the bacterial orders with highest relative abundances in the samples in order to choose which orders to display.

## There are different ways to do this:

### 1) One can look at the top20 bacterial orders with highest mean relative abundance in the different soil-compartment combinations for Lotus and Hordeum separately.
##### With this option one would however end up with a large number of different bacterial orders.
##### (Even when filtering for orders with e.g. mean relative abundance of min. 1%).

### 2) One can go by the top20 bacterial orders in either the rhizosphere OR the root of Lotus or Hordeum (no separation by soil types, only by compartment).
##### So top20 orders in rhizosphere and top20 in root for Lotus or Hordeum would be identified separately.
##### The unique orders would be combined, and one could additionally filter for RA of min. 1% in at least one of the compartments in one of the plants.

## NOTE: The Lotus nodule compartment was not taken into account, as it is mainly dominated by Rhizobiales and top20 would represent orders with very small RA.

# Step 1: Filter for rhizosphere and root compartments only.
Lotus_RR <- Lotus_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

Hordeum_RR <- Hordeum_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Step 2: Calculate the mean RA per bacterial order across all samples of each compartment (soil samples combined).
Lotus_meanRA <- Lotus_RR %>%
  group_by(Compartment, Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

Hordeum_meanRA <- Hordeum_RR %>%
  group_by(Compartment, Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

# Step 3: Identify the top 20 orders in each compartment by the mean RA.
Lotus_top20_orders_per_comp <- Lotus_meanRA %>%
  group_by(Compartment) %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20) %>%
  ungroup()

Hordeum_top20_orders_per_comp <- Hordeum_meanRA %>%
  group_by(Compartment) %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20) %>%
  ungroup()

# Step 4: Make a unique list of orders across the two compartments.
Lotus_unique_top20_orders <- Lotus_top20_orders_per_comp %>%
  distinct(Order) %>%
  pull(Order)

Hordeum_unique_top20_orders <- Hordeum_top20_orders_per_comp %>%
  distinct(Order) %>%
  pull(Order)

# Step 5: Check which of these orders have a mean RA >= 0.01 in either rhizosphere or root.
Lotus_top_orders_above1perc <- Lotus_meanRA %>%
  filter(Order %in% Lotus_unique_top20_orders & MeanRA >= 0.01) %>%
  distinct(Order) %>%
  pull(Order)

Hordeum_top_orders_above1perc <- Lotus_meanRA %>%
  filter(Order %in% Hordeum_unique_top20_orders & MeanRA >= 0.01) %>%
  distinct(Order) %>%
  pull(Order)

# Step 6: Check how many unique orders these would be to display when combining
# Lotus and Hordeum top20 rhizosphere and/or root with >=1% RA.
combined_top_orders_1 <- unique(c(Lotus_top_orders_above1perc, Hordeum_top_orders_above1perc))
length(combined_top_orders_1)

### 3) One can go by the top20 bacterial orders across all samples, identified for 
### Lotus and Hordeum separately.

##### One would not look for top20 separately for soil or compartment, but average
##### it across all samples. Then one could filter for orders with mean RA of min.
##### 1%, and combine the unique orders for Lotus and Hordeum and display those.

## NOTE: Lotus nodule compartment again not taken into account, so average in 
## rhizosphere and root calculated. I did test taking the nodule compartment 
## into account, but the result looked exactly the same, apart from loosing
## one group when nodule samples were included. Therefore I kept it excluded.

# Step 1: Keep only the rhizosphere and root compartments.
Lotus_filtered <- Lotus_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

Hordeum_filtered <- Hordeum_order_summary %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Step 2: Calculate the overall mean RA per order across these samples.
Lotus_meanRA <- Lotus_filtered %>%
  group_by(Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

Hordeum_meanRA <- Hordeum_filtered %>%
  group_by(Order) %>%
  summarise(MeanRA = mean(RA, na.rm = TRUE), .groups = "drop")

# Step 3: Identify the top 20 orders by overall mean RA.
Lotus_top20_orders <- Lotus_meanRA %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20)

Hordeum_top20_orders <- Hordeum_meanRA %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 20)

# Step 4: Keep only those with mean RA ≥ 0.01.
Lotus_top20_above1perc <- Lotus_top20_orders %>%
  filter(MeanRA >= 0.01) %>%
  pull(Order)

Hordeum_top20_above1perc <- Hordeum_top20_orders %>%
  filter(MeanRA >= 0.01) %>%
  pull(Order)

# Step 5: Combine the orders for Lotus and Hordeum and show how many unique orders 
# would be displayed.
combined_top_orders_2 <- unique(c(Lotus_top20_above1perc, Hordeum_top20_above1perc))
length(combined_top_orders_2)

####### I WILL GO AHEAD WITH OPTION 3 (combined_top_orders_2).

# We will work with the selected bacterial orders, go back to the original 
# dataframes and rename all other bacterial orders "Other".
Lotus_df.long <- Lotus_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders_2, Order, "Other"))

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Order = if_else(Order %in% combined_top_orders_2, Order, "Other"))

# We then add a new column "Plant" to both dataframes and combine them.
Lotus_df.long <- Lotus_df.long %>%
  mutate(Plant = "Lotus")

Hordeum_df.long <- Hordeum_df.long %>%
  mutate(Plant = "Hordeum")

combined_df <- bind_rows(Lotus_df.long, Hordeum_df.long)

# We now summarise the mean RA per plant-compartment-soil combination.
df.sample_order <- combined_df %>%
  group_by(sampleID, Plant, Compartment, Soil, Order) %>%
  summarise(RA=sum(RA), .groups="drop")

df.mean_order <- df.sample_order %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  summarise(RA=mean(RA), .groups="drop") %>%
  mutate(Order=factor(Order, levels=c(sort(unique(Order[Order!="Other"])), "Other")),
         Plant=factor(Plant, levels=c("Lotus","Hordeum")),
         Compartment=factor(Compartment, levels=c("Rhizosphere","Root","Nodules")),
         Soil=factor(Soil, levels=c("NPK","PK","UF")))


# Set the main theme for the stacked barplot.
main_theme <- theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  panel.border=element_rect(colour="black", fill=NA, linewidth=1),
  axis.line.x=element_line(color="black"),
  axis.line.y=element_line(color="black"),
  axis.ticks=element_line(color="black"),
  axis.text=element_text(size = 6, color="black"),
  legend.text=element_text(size = 6, color="black"),
  legend.key=element_blank(),
  axis.title.y=element_text(size = 6),
  legend.position="right",
  legend.background=element_blank(),
  text=element_text(family="sans", size = 6, color="black")
)

# Set the colours for the bacterial orders.
colors <- read.table("../../../0_files/Bacterial_order_colors.csv", header = T, sep = ",", comment.char = "")

# Make the stacked barplot.
p1 <- ggplot(df.mean_order, aes(x=Soil, y=RA, fill=Order)) +
  geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values = colors$Color, breaks = colors$Order) +
  scale_y_continuous(expand=c(0,0)) +
  main_theme +
  ylab("Mean relative abundance") +
  labs(fill="Bacterial order") +
  xlab("") +
  guides(fill=guide_legend(nrow=21)) +
  facet_nested(~ Plant + Compartment, scales="free_x", space="free_x") +
  theme(
    axis.text.x = element_text(size = 6, color="black", angle=0, vjust=1),
    strip.text.x = element_text(size = 6, face="bold"),
    legend.key.size = unit(0.25, 'cm'),
    legend.margin = margin(l = -8)
  )

p1

# Save the plot.
ggsave("LotusHordeum_Askov_WT_stackedbp_top20_meanRA.pdf", p1, width = 10, height = 6, unit = "cm")
saveRDS(p1, file = "LotusHordeum_Askov_WT_stackedbp_top20_meanRA.rds")
saveRDS(p1, file = "../7_final_figures/LotusHordeum_Askov_WT_stackedbp_top20_meanRA.rds")

# Make a heatmap that displays these data and includes info on significance 
# (differences between soil type in each plant-compartment combination.)
# Prepare data for the heatmap by removing nodule samples.
df.heatmap <- df.sample_order %>%
  filter(Compartment %in% c("Rhizosphere", "Root"))

# Calculate the mean RA per plant-compartment-soil-order.
df.heatmap <- df.heatmap %>%
  group_by(Plant, Compartment, Soil, Order) %>%
  summarise(RA = mean(RA, na.rm = TRUE), .groups="drop")

# Set factor levels and arrange the orders alphabetical, with "Other" last.
df.heatmap <- df.heatmap %>%
  mutate(Order = factor(Order,
                        levels = c(sort(unique(Order[Order != "Other"])), "Other")),
         Plant = factor(Plant, levels=c("Lotus","Hordeum")),
         Compartment = factor(Compartment, levels=c("Rhizosphere","Root")),
         Soil = factor(Soil, levels=c("NPK","PK","UF")))

# Perform significance analysis.
## Initialise a list for the significance letters.
letters_list <- list()
plants <- unique(df.heatmap$Plant)
comps  <- unique(df.heatmap$Compartment)
ords  <- unique(df.heatmap$Order)
cond.grid <- expand.grid(Plant = plants, Compartment = comps, Order = ords)
j <- 1

## Loop over each plant-compartment-order combination.
for(i in 1:nrow(cond.grid)){
  df_samples <- df.sample_order %>% 
        filter(
          Plant == cond.grid$Plant[i], 
          Compartment == cond.grid$Compartment[i],
          Order == cond.grid$Order[i]
        )
  df_samples$Soil <- factor(df_samples$Soil, levels = c("UF", "NPK", "PK"))
  a <- aov(RA ~ Soil, data = df_samples)
  anv <- anova(a)
  p_val <- anv$`Pr(>F)`[1]
  if(p_val < 0.05){
    tk <- TukeyHSD(a)
    mcletters <- multcompLetters(tk$Soil[,"p adj"])
    letters <- mcletters$Letters
    df_res <- data.frame(
      Plant = cond.grid$Plant[i],
      Compartment = cond.grid$Compartment[i],
      Order = cond.grid$Order[i],
      Soil = names(letters),
      letter = letters,
      stringsAsFactors = FALSE
    )
    letters_list[[j]] <- df_res
    j <- j+1
  }
}

# ## Loop over each plant-compartment-order.
# for(pl in plants){
#   for(comp in comps){
#     df_sub <- df.heatmap %>% filter(Plant==pl, Compartment==comp)
#     for(ord in unique(df_sub$Order)){
#       df_ord <- df_sub %>% filter(Order==ord)
#       df_samples <- df.sample_order %>% 
#         filter(Plant==pl, Compartment==comp, Order==ord)
#       soil_counts <- table(df_samples$Soil)
#       # Only test if at least 2 soils have >1 sample in original df.sample_order.
#       if(length(soil_counts) > 1 & all(soil_counts >= 2)){
#         kw <- kruskal.test(RA ~ Soil, data=df_samples)
#         if(!is.na(kw$p.value) & kw$p.value < 0.05){
#           # Perform Dunn post-hoc.
#           dunn_res <- dunnTest(RA ~ Soil, data=df_samples, method="bh")$res
#           # Convert p-adj to letters.
#           pvals_named <- setNames(dunn_res$P.adj, gsub(" ", "", dunn_res$Comparison))
#           cld <- multcompLetters(pvals_named, threshold=0.05)
#           # Map letters to soil levels.
#           soil_levels <- levels(df.heatmap$Soil)
#           letter_vec <- rep(NA, length(soil_levels))
#           names(letter_vec) <- soil_levels
#           letter_vec[names(cld$Letters)] <- cld$Letters
          
#           letters_list[[paste(pl, comp, ord, sep="_")]] <- data.frame(
#             Plant = pl,
#             Compartment = comp,
#             Order = ord,
#             Soil = soil_levels,
#             letter = letter_vec,
#             stringsAsFactors = FALSE
#           )
#         }
#       }
#     }
#   }
# }

## Combine all letters into a dataframe.
df_letters <- bind_rows(letters_list)

# Join the significance letters with the heatmap data.
df.plot <- df.heatmap %>%
  left_join(df_letters, by=c("Plant","Compartment","Order","Soil"))
# df.plot <- df.heatmap %>%
#   left_join(df_letters, by=c("Plant","Compartment","Order","Soil")) %>%
#   mutate(
#     Plant = factor(Plant, levels=c("Lotus","Hordeum")),
#     Compartment = factor(Compartment, levels=c("Rhizosphere","Root")),
#     Order = factor(Order, levels = c(sort(unique(Order[Order != "Other"])), "Other"))
#   )

# Define the axis breaks and colours.
breaks <- c(0, 0.005, 0.052, 0.052001, 0.15999, 0.16, 0.34, 0.64)
colors <- c("#1F78B4", "#A6CEE3", "white","#FFFF99",
            "#FF7F00", "#FB9A99", "#E31A1C", "#902121")

# Rescale breaks to 0-1 for gradient.
values <- rescale(breaks, to = c(0,1))

letter_keep <- df.plot %>% 
  group_by(Plant, Compartment, Order) %>% 
  summarise(remove = all(letter == "a")) %>%
  ungroup()

letter_keep$remove[is.na(letter_keep$remove)] <- F
Opt <- expand.grid(Plant = c("Lotus", "Hordeum"), 
                   Compartment = c("Root", "Rhizosphere"))

for(i in 1:nrow(Opt)){
  letter_keep_sub <- letter_keep[
    letter_keep$Plant == Opt$Plant[i] & 
    letter_keep$Compartment == Opt$Compartment[i],
  ]
  orders_remove_letters <- as.character(letter_keep_sub$Order[letter_keep_sub$remove])
  df.plot$letter[
    df.plot$Plant == Opt$Plant[i] & 
    df.plot$Compartment == Opt$Compartment[i] &
    df.plot$Order %in% orders_remove_letters
  ] <- NA
}

df.plot <- df.plot %>% 
  filter(!(Order %in% c("Other", "Unknown"))) %>% 
  mutate(Order = droplevels(Order))

# Then make the final heatmap.
p_heatmap <- ggplot(df.plot, aes(x=Soil, y=Order, fill=RA)) +
  geom_tile(color="grey50") +
  geom_text(aes(label=letter), na.rm=TRUE, size = 6/.pt) +
  scale_fill_gradientn(
    colors = colors,
    values = values,
    limits = c(0, max(df.plot$RA, na.rm=TRUE)),
    name = "Relative abundance"
  ) +
  scale_y_discrete(
    limits = rev(levels(df.plot$Order)),
    position = "right"
  ) +
  guides(fill = guide_colorbar(
    title.position = "right",
    barwidth = 12.5,
    barheight = 1
  ))+
  main_theme +
  facet_nested(~ Plant + Compartment, scales="free_x", space="free_x") +
  xlab(NULL) +
  ylab("Bacterial order") +
  theme(
    axis.text.x = element_text(size = 6, angle=0, vjust=1, hjust=0.5, 
                               colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 6, colour = "black"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill="grey90", color=NA),
    strip.text = element_text(size = 6, face="bold"),
    legend.text = element_text(size = 6, colour = "black"),
    legend.title = element_text(size = 6, colour = "black"),
    legend.position = "bottom",
    plot.margin = margin(r = 10, l = 20)
  )

p_heatmap

# Save the heatmap.
ggsave("LotusHordeum_Askov_WT_orders_heatmap.pdf", p_heatmap, width = 12, height = 6, unit = "cm")
saveRDS(p_heatmap, file = "LotusHordeum_Askov_WT_orders_heatmap.rds")
saveRDS(p_heatmap, file = "../7_final_figures/LotusHordeum_Askov_WT_orders_heatmap.rds")
