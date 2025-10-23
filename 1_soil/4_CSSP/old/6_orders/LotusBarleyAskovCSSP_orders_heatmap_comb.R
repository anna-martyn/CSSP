# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load files.
Lotus_design.file <- paste("Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", sep="")
Lotus_taxonomy.file <- paste("LotusSep_exclUFnew_10_4_silva138_taxonomy-3.txt", sep="")
Lotus_otu_table.file <- paste("LotusSep_exclUFnew_ASVtable_10_4_nospike.txt", sep="")

Lotus_design <- read.table(Lotus_design.file, header=T, sep="\t")
Lotus_otu_table <- read.table(Lotus_otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
Lotus_taxonomy <- read.table(Lotus_taxonomy.file, sep="\t", header=T, fill=T)

Barley_design.file <- paste("BarleyCSSP_Askov_reseq_metadata-2.txt", sep="")
Barley_taxonomy.file <- paste("Barley_Askov_Rep_10_4_taxonomy.txt", sep="")
Barley_otu_table.file <- paste("BarleyCSSP_Askov_reseq_ASVtable_10_4.txt", sep="")

Barley_design <- read.table(Barley_design.file, header=T, sep="\t")
Barley_otu_table <- read.table(Barley_otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
Barley_taxonomy <- read.table(Barley_taxonomy.file, sep="\t", header=T, fill=T)


# Load required packages.
library(dplyr)
library(data.table)
library(ggplot2)
library(multcompView)
library(openxlsx)

# Clean-up taxonomy file layout.
Lotus_taxa <- lapply(as.list(Lotus_taxonomy$Taxon),
               function(x) x %>% strsplit(split = "; ") %>% unlist() %>% substr(start = 4, stop = 100))

Lotus_taxa <- lapply(Lotus_taxa, function(x) c( x, rep("Unknown", 7 - length(x)) ))
Lotus_taxa <- lapply(Lotus_taxa, function(x) 
  data.table(Kingdom = x[1],
             Phylum = x[2],
             Class = x[3],
             Order = x[4],
             Family = x[5],
             Genus = x[6],
             Species = x[7])
)
Lotus_taxa <- rbindlist(Lotus_taxa)
Lotus_taxonomy <- data.table(ASV = Lotus_taxonomy$`Feature.ID`, Lotus_taxa, Confidence = Lotus_taxonomy$Confidence)

Barley_taxa <- lapply(as.list(Barley_taxonomy$Taxon),
                     function(x) x %>% strsplit(split = "; ") %>% unlist() %>% substr(start = 4, stop = 100))

Barley_taxa <- lapply(Barley_taxa, function(x) c( x, rep("Unknown", 7 - length(x)) ))
Barley_taxa <- lapply(Barley_taxa, function(x) 
  data.table(Kingdom = x[1],
             Phylum = x[2],
             Class = x[3],
             Order = x[4],
             Family = x[5],
             Genus = x[6],
             Species = x[7])
)
Barley_taxa <- rbindlist(Barley_taxa)
Barley_taxonomy <- data.table(ASV = Barley_taxonomy$`ASVid`, Barley_taxa, Confidence = Barley_taxonomy$Confidence)

# re-order data matrices
idx <- Lotus_design$SampleID %in% colnames(Lotus_otu_table)
Lotus_design <- Lotus_design[idx, ]

idx <- match(Lotus_design$SampleID, colnames(Lotus_otu_table))
Lotus_otu_table <- Lotus_otu_table[, idx]

idx <- rownames(Lotus_otu_table) %in% Lotus_taxonomy$ASV
Lotus_otu_table <- Lotus_otu_table[idx, ]

idx <- match(rownames(Lotus_otu_table), Lotus_taxonomy$ASV)
Lotus_taxonomy <- Lotus_taxonomy[idx, ]

idx <- Barley_design$Sample_ID %in% colnames(Barley_otu_table)
Barley_design <- Barley_design[idx, ]

idx <- match(Barley_design$Sample_ID, colnames(Barley_otu_table))
Barley_otu_table <- Barley_otu_table[, idx]

idx <- rownames(Barley_otu_table) %in% Barley_taxonomy$ASV
Barley_otu_table <- Barley_otu_table[idx, ]

idx <- match(rownames(Barley_otu_table), Barley_taxonomy$ASV)
Barley_taxonomy <- Barley_taxonomy[idx, ]

# otu_table normalization
Lotus_otu_table_norm <- apply(Lotus_otu_table, 2, function(x) x / sum(x))
Lotus_df <- as.data.frame(Lotus_otu_table_norm)
Lotus_df$ASV <- row.names(Lotus_df)
row.names(Lotus_df) <- NULL

Barley_otu_table_norm <- apply(Barley_otu_table, 2, function(x) x / sum(x))
Barley_df <- as.data.frame(Barley_otu_table_norm)
Barley_df$ASV <- row.names(Barley_df)
row.names(Barley_df) <- NULL

### reshape data: the aim of this step is to transfer wide data to long data(check online the difference of wide and long dataset in r)
Lotus_df.long <- reshape(Lotus_df, varying = 1:204, direction = "long", idvar = 'ASV', timevar="sampleID", v.names=c("RA"), 
                   times=c(colnames(Lotus_df[,1:204])), sep="")
row.names(Lotus_df.long) <- NULL
Lotus_df.long <- as.data.frame(Lotus_df.long)

Barley_df.long <- reshape(Barley_df, varying = 1:191, direction = "long", idvar = 'ASV', timevar="sampleID", v.names=c("RA"), 
                         times=c(colnames(Barley_df[,1:191])), sep="")
row.names(Barley_df.long) <- NULL
Barley_df.long <- as.data.frame(Barley_df.long)

# Add taxonomy information (order).
Lotus_df.long_order <- merge(Lotus_df.long, Lotus_taxonomy[, c("ASV", "Order")], by = "ASV", all.x = TRUE)

Barley_df.long_order <- merge(Barley_df.long, Barley_taxonomy[, c("ASV", "Order")], by = "ASV", all.x = TRUE)

###  make a dataframe that sum RA by tax and sampleID
Lotus_df.long_order <- Lotus_df.long_order %>% select (2:4) ### select the column 2,3,4 from df.long dataframe

Lotus_df.long_order <- Lotus_df.long_order %>%
  group_by_(.dots = c("Order", "sampleID")) %>%   ### this is to summarize RA by both tax and sampleID
  summarise_all(funs(sum))

Barley_df.long_order <- Barley_df.long_order %>% select (2:4) ### select the column 2,3,4 from df.long dataframe

Barley_df.long_order <- Barley_df.long_order %>%
  group_by_(.dots = c("Order", "sampleID")) %>%   ### this is to summarize RA by both tax and sampleID
  summarise_all(funs(sum))

### Add genotype information for df.long.
Lotus_df3 <- data.frame(sampleID = Lotus_design$SampleID, Genotype=Lotus_design$Genotype, Compartment = Lotus_design$Compartment, Plant=Lotus_design$Plant, Soil=Lotus_design$Soil)

Lotus_df.long_order <- full_join(Lotus_df.long_order, Lotus_df3, by="sampleID")


idx <- Lotus_df.long_order$Genotype%in% c("WT","symrk","ccamk","nsp1","nsp2")
Lotus_df.long_order <- Lotus_df.long_order[idx,]

idx <- Lotus_df.long_order$Soil%in% c("NPK")
Lotus_df.long_order_NPK <- Lotus_df.long_order[idx,]

idx <- Lotus_df.long_order_NPK$Compartment%in% c("rhizosphere")
Lotus_df.long_order_NPK_rhizo <- Lotus_df.long_order_NPK[idx,]

idx <- Lotus_df.long_order_NPK$Compartment%in% c("root")
Lotus_df.long_order_NPK_root <- Lotus_df.long_order_NPK[idx,]


idx <- Lotus_df.long_order$Soil%in% c("PK")
Lotus_df.long_order_PK <- Lotus_df.long_order[idx,]

idx <- Lotus_df.long_order_PK$Compartment%in% c("rhizosphere")
Lotus_df.long_order_PK_rhizo <- Lotus_df.long_order_PK[idx,]

idx <- Lotus_df.long_order_PK$Compartment%in% c("root")
Lotus_df.long_order_PK_root <- Lotus_df.long_order_PK[idx,]

idx <- Lotus_df.long_order$Soil%in% c("UF")
Lotus_df.long_order_UF <- Lotus_df.long_order[idx,]

idx <- Lotus_df.long_order_UF$Compartment%in% c("rhizosphere")
Lotus_df.long_order_UF_rhizo <- Lotus_df.long_order_UF[idx,]

idx <- Lotus_df.long_order_UF$Compartment%in% c("root")
Lotus_df.long_order_UF_root <- Lotus_df.long_order_UF[idx,]


Barley_df3 <- data.frame(sampleID = Barley_design$Sample_ID, Genotype=Barley_design$Genotype, Compartment = Barley_design$Compartment, Plant=Barley_design$Plant, Soil=Barley_design$Soil)

Barley_df.long_order <- full_join(Barley_df.long_order, Barley_df3, by="sampleID")


idx <- Barley_df.long_order$Genotype%in% c("WT","symrk","ccamk","nsp1","nsp2")
Barley_df.long_order <- Barley_df.long_order[idx,]

idx <- Barley_df.long_order$Soil%in% c("NPK")
Barley_df.long_order_NPK <- Barley_df.long_order[idx,]

idx <- Barley_df.long_order_NPK$Compartment%in% c("rhizosphere")
Barley_df.long_order_NPK_rhizo <- Barley_df.long_order_NPK[idx,]

idx <- Barley_df.long_order_NPK$Compartment%in% c("root")
Barley_df.long_order_NPK_root <- Barley_df.long_order_NPK[idx,]


idx <- Barley_df.long_order$Soil%in% c("PK")
Barley_df.long_order_PK <- Barley_df.long_order[idx,]

idx <- Barley_df.long_order_PK$Compartment%in% c("rhizosphere")
Barley_df.long_order_PK_rhizo <- Barley_df.long_order_PK[idx,]

idx <- Barley_df.long_order_PK$Compartment%in% c("root")
Barley_df.long_order_PK_root <- Barley_df.long_order_PK[idx,]

idx <- Barley_df.long_order$Soil%in% c("UF")
Barley_df.long_order_UF <- Barley_df.long_order[idx,]

idx <- Barley_df.long_order_UF$Compartment%in% c("rhizosphere")
Barley_df.long_order_UF_rhizo <- Barley_df.long_order_UF[idx,]

idx <- Barley_df.long_order_UF$Compartment%in% c("root")
Barley_df.long_order_UF_root <- Barley_df.long_order_UF[idx,]

# Caluclate mean relative abundance of all orders taking all genotypes into account in the individual soil-fraction combinations.
# Then check which ones have a mean relative abundance >1% in any of the conditions.

## Calculate the mean RA.
Lotus_process_df <- function(Lotus_df) {
  Lotus_df %>%
    group_by(Order) %>%
    summarize(Mean = mean(RA, na.rm = TRUE)) %>%
    arrange(desc(Mean)) %>%
    filter(Mean >= 0.01) %>%
    as.data.frame()
}

Barley_process_df <- function(Barley_df) {
  Barley_df %>%
    group_by(Order) %>%
    summarize(Mean = mean(RA, na.rm = TRUE)) %>%
    arrange(desc(Mean)) %>%
    filter(Mean >= 0.01) %>%
    as.data.frame()
}

# Apply the function to each dataframe
Lotus_mean_RA_order_NPK_rhizo_0.01 <- Lotus_process_df(Lotus_df.long_order_NPK_rhizo)
Lotus_mean_RA_order_PK_rhizo_0.01  <- Lotus_process_df(Lotus_df.long_order_PK_rhizo)
Lotus_mean_RA_order_UF_rhizo_0.01  <- Lotus_process_df(Lotus_df.long_order_UF_rhizo)
Lotus_mean_RA_order_NPK_root_0.01  <- Lotus_process_df(Lotus_df.long_order_NPK_root)
Lotus_mean_RA_order_PK_root_0.01   <- Lotus_process_df(Lotus_df.long_order_PK_root)
Lotus_mean_RA_order_UF_root_0.01   <- Lotus_process_df(Lotus_df.long_order_UF_root)

Barley_mean_RA_order_NPK_rhizo_0.01 <- Barley_process_df(Barley_df.long_order_NPK_rhizo)
Barley_mean_RA_order_PK_rhizo_0.01  <- Barley_process_df(Barley_df.long_order_PK_rhizo)
Barley_mean_RA_order_UF_rhizo_0.01  <- Barley_process_df(Barley_df.long_order_UF_rhizo)
Barley_mean_RA_order_NPK_root_0.01  <- Barley_process_df(Barley_df.long_order_NPK_root)
Barley_mean_RA_order_PK_root_0.01   <- Barley_process_df(Barley_df.long_order_PK_root)
Barley_mean_RA_order_UF_root_0.01   <- Barley_process_df(Barley_df.long_order_UF_root)

## Check all unique orders above 1% mean relative abundance in any of the created dataframes.
Lotus_unique_orders <- bind_rows(Lotus_mean_RA_order_NPK_rhizo_0.01, Lotus_mean_RA_order_PK_rhizo_0.01, Lotus_mean_RA_order_UF_rhizo_0.01, Lotus_mean_RA_order_NPK_root_0.01, Lotus_mean_RA_order_PK_root_0.01, Lotus_mean_RA_order_UF_root_0.01) %>%
  distinct(Order) %>%
  pull(Order)

Barley_unique_orders <- bind_rows(Barley_mean_RA_order_NPK_rhizo_0.01, Barley_mean_RA_order_PK_rhizo_0.01, Barley_mean_RA_order_UF_rhizo_0.01, Barley_mean_RA_order_NPK_root_0.01, Barley_mean_RA_order_PK_root_0.01, Barley_mean_RA_order_UF_root_0.01) %>%
  distinct(Order) %>%
  pull(Order)

both_unique_orders <- bind_rows(Lotus_mean_RA_order_NPK_rhizo_0.01, Lotus_mean_RA_order_PK_rhizo_0.01, Lotus_mean_RA_order_UF_rhizo_0.01, Lotus_mean_RA_order_NPK_root_0.01, Lotus_mean_RA_order_PK_root_0.01, Lotus_mean_RA_order_UF_root_0.01,Barley_mean_RA_order_NPK_rhizo_0.01, Barley_mean_RA_order_PK_rhizo_0.01, Barley_mean_RA_order_UF_rhizo_0.01, Barley_mean_RA_order_NPK_root_0.01, Barley_mean_RA_order_PK_root_0.01, Barley_mean_RA_order_UF_root_0.01) %>%
  distinct(Order) %>%
  pull(Order)
  
# Now filter the original dataframe to only keep the orders of interest (RA >1%),a nd only focus on rhizosphere and root fractions.
Lotus_df_filtered <- Lotus_df.long_order[Lotus_df.long_order$Order %in% both_unique_orders, ]

Lotus_df_filtered <- Lotus_df_filtered %>%
  filter(Order != "Unknown")

Lotus_df_filtered <- Lotus_df_filtered[Lotus_df_filtered$Compartment %in% c("rhizosphere", "root"), ]

Lotus_df_filtered$Genotype <- factor(Lotus_df_filtered$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))

Barley_df_filtered <- Barley_df.long_order[Barley_df.long_order$Order %in% both_unique_orders, ]

Barley_df_filtered <- Barley_df_filtered %>%
  filter(Order != "Unknown")

Barley_df_filtered <- Barley_df_filtered[Barley_df_filtered$Compartment %in% c("rhizosphere", "root"), ]

Barley_df_filtered$Genotype <- factor(Barley_df_filtered$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))

# df_filtered$Order <- factor(df_filtered$Order, levels = sort(unique(df_filtered$Order)))
Lotus_df_filtered$Order <- factor(Lotus_df_filtered$Order, levels = rev(sort(unique(Lotus_df_filtered$Order))))

Barley_df_filtered$Order <- factor(Barley_df_filtered$Order, levels = rev(sort(unique(Barley_df_filtered$Order))))

# Now make plot.
## Define breaks and colours for RA legend:
breaks <- c(seq(0.00, 0.05, by = 0.01), seq(0.06, 0.20, by = 0.04), seq(0.21, 1.0, by = 0.1))
colors <- c("darkblue", "deepskyblue", "white", "#feedb4", "#F09D00", "#FFB3B2", "#ba0319", "#902121")

# Plot the heatmap
main_theme <- theme(axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="right",
                    # legend.background=element_blank(),
                    # legend.key=element_blank(),
                    text=element_text(family="sans"))

p1 <- ggplot(Lotus_df_filtered, aes(x = Genotype, y = Order, fill = RA)) +
  geom_tile() +
  scale_fill_gradientn(colors = colors, values = scales::rescale(breaks)) +
  facet_grid(. ~ Compartment + Soil) +
  labs(x = "Genotype", y = "Order", fill = "RA") +
  main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, hjust=1))

ggsave(paste("LotusCSSP_heatmap_orders_above_equal1perc_comb.png", sep=""), p1, width=14, height=12)
ggsave(paste("LotusCSSP_heatmap_orders_above_equal1perc_comb.pdf", sep=""), p1, width=14, height=12)

p2 <- ggplot(Barley_df_filtered, aes(x = Genotype, y = Order, fill = RA)) +
  geom_tile() +
  scale_fill_gradientn(colors = colors, values = scales::rescale(breaks)) +
  facet_grid(. ~ Compartment + Soil) +
  labs(x = "Genotype", y = "Order", fill = "RA") +
  main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 45, hjust=1))

ggsave(paste("BarleyCSSP_heatmap_orders_above_equal1perc_comb.png", sep=""), p2, width=14, height=12)
ggsave(paste("BarleyCSSP_heatmap_orders_above_equal1perc_comb.pdf", sep=""), p2, width=14, height=12)

# ##########SIGNIFICANCE TESTING##################

####Lotus 
Lotus_df_filtered_rhizo <- Lotus_df_filtered[Lotus_df_filtered$Compartment %in% c("rhizosphere"), ]
Lotus_df_filtered_root <- Lotus_df_filtered[Lotus_df_filtered$Compartment %in% c("root"), ]

Lotus_df_filtered_rhizo_NPK <- Lotus_df_filtered_rhizo[Lotus_df_filtered_rhizo$Soil %in% c("NPK"), ]
Lotus_df_filtered_rhizo_PK <- Lotus_df_filtered_rhizo[Lotus_df_filtered_rhizo$Soil %in% c("PK"), ]
Lotus_df_filtered_rhizo_UF <- Lotus_df_filtered_rhizo[Lotus_df_filtered_rhizo$Soil %in% c("UF"), ]

Lotus_df_filtered_root_NPK <- Lotus_df_filtered_root[Lotus_df_filtered_root$Soil %in% c("NPK"), ]
Lotus_df_filtered_root_PK <- Lotus_df_filtered_root[Lotus_df_filtered_root$Soil %in% c("PK"), ]
Lotus_df_filtered_root_UF <- Lotus_df_filtered_root[Lotus_df_filtered_root$Soil %in% c("UF"), ]


##################NPK_rhizo
orders <- unique(Lotus_df_filtered_rhizo_NPK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Lotus_df_filtered_rhizo_NPK$Order)) {
  idx <- Lotus_df_filtered_rhizo_NPK$Order %in% c(order)
  Lotus_df_filtered <- Lotus_df_filtered_rhizo_NPK[idx, ]
  if (nrow(Lotus_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Lotus_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_rhizo_NPK_ANOVA_comb.xlsx", row.names = FALSE)

##################PK_rhizo
orders <- unique(Lotus_df_filtered_rhizo_PK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Lotus_df_filtered_rhizo_PK$Order)) {
  idx <- Lotus_df_filtered_rhizo_PK$Order %in% c(order)
  Lotus_df_filtered <- Lotus_df_filtered_rhizo_PK[idx, ]
  if (nrow(Lotus_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Lotus_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_rhizo_PK_ANOVA_comb.xlsx", row.names = FALSE)

##################UF_rhizo
orders <- unique(Lotus_df_filtered_rhizo_UF$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Lotus_df_filtered_rhizo_UF$Order)) {
  idx <- Lotus_df_filtered_rhizo_UF$Order %in% c(order)
  Lotus_df_filtered <- Lotus_df_filtered_rhizo_UF[idx, ]
  if (nrow(Lotus_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Lotus_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_rhizo_UF_ANOVA_comb.xlsx", row.names = FALSE)

##################NPK_root
orders <- unique(Lotus_df_filtered_root_NPK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Lotus_df_filtered_root_NPK$Order)) {
  idx <- Lotus_df_filtered_root_NPK$Order %in% c(order)
  Lotus_df_filtered <- Lotus_df_filtered_root_NPK[idx, ]
  if (nrow(Lotus_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Lotus_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_root_NPK_ANOVA_comb.xlsx", row.names = FALSE)

##################PK_root
orders <- unique(Lotus_df_filtered_root_PK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Lotus_df_filtered_root_PK$Order)) {
  idx <- Lotus_df_filtered_root_PK$Order %in% c(order)
  Lotus_df_filtered <- Lotus_df_filtered_root_PK[idx, ]
  if (nrow(Lotus_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Lotus_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_root_PK_ANOVA_comb.xlsx", row.names = FALSE)

##################UF_root
orders <- unique(Lotus_df_filtered_root_UF$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable) {
  if (is.null(pairwise[[variable]])) {
    return(data.frame(Genotype = genotype_order, Letters = NA))  # Return NAs if no results
  }

  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  if (any(is.na(Tukey.levels))) {
    return(data.frame(Genotype = genotype_order, Letters = NA))  # Handle missing values
  }

  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

# Define genotype order
genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")

# Initialize final dataframe
final_df <- data.frame(Genotype = genotype_order)

for (order in orders) {
  idx <- Lotus_df_filtered_root_UF$Order %in% c(order)
  Lotus_df_filtered <- Lotus_df_filtered_root_UF[idx, ]
  if (nrow(Lotus_df_filtered) < length(genotype_order)) {
    next
  }
  ano <- aov(RA ~ Genotype, data = Lotus_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_root_UF_ANOVA_comb.xlsx", row.names = FALSE)












####Barley 
Barley_df_filtered_rhizo <- Barley_df_filtered[Barley_df_filtered$Compartment %in% c("rhizosphere"), ]
Barley_df_filtered_root <- Barley_df_filtered[Barley_df_filtered$Compartment %in% c("root"), ]

Barley_df_filtered_rhizo_NPK <- Barley_df_filtered_rhizo[Barley_df_filtered_rhizo$Soil %in% c("NPK"), ]
Barley_df_filtered_rhizo_PK <- Barley_df_filtered_rhizo[Barley_df_filtered_rhizo$Soil %in% c("PK"), ]
Barley_df_filtered_rhizo_UF <- Barley_df_filtered_rhizo[Barley_df_filtered_rhizo$Soil %in% c("UF"), ]

Barley_df_filtered_root_NPK <- Barley_df_filtered_root[Barley_df_filtered_root$Soil %in% c("NPK"), ]
Barley_df_filtered_root_PK <- Barley_df_filtered_root[Barley_df_filtered_root$Soil %in% c("PK"), ]
Barley_df_filtered_root_UF <- Barley_df_filtered_root[Barley_df_filtered_root$Soil %in% c("UF"), ]


##################NPK_rhizo
orders <- unique(Barley_df_filtered_rhizo_NPK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Barley_df_filtered_rhizo_NPK$Order)) {
  idx <- Barley_df_filtered_rhizo_NPK$Order %in% c(order)
  Barley_df_filtered <- Barley_df_filtered_rhizo_NPK[idx, ]
  if (nrow(Barley_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Barley_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "BarleyCSSP_rhizo_NPK_ANOVA_comb.xlsx", row.names = FALSE)

##################PK_rhizo
orders <- unique(Barley_df_filtered_rhizo_PK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Barley_df_filtered_rhizo_PK$Order)) {
  idx <- Barley_df_filtered_rhizo_PK$Order %in% c(order)
  Barley_df_filtered <- Barley_df_filtered_rhizo_PK[idx, ]
  if (nrow(Barley_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Barley_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "BarleyCSSP_rhizo_PK_ANOVA_comb.xlsx", row.names = FALSE)

##################UF_rhizo
orders <- unique(Barley_df_filtered_rhizo_UF$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Barley_df_filtered_rhizo_UF$Order)) {
  idx <- Barley_df_filtered_rhizo_UF$Order %in% c(order)
  Barley_df_filtered <- Barley_df_filtered_rhizo_UF[idx, ]
  if (nrow(Barley_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Barley_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "BarleyCSSP_rhizo_UF_ANOVA_comb.xlsx", row.names = FALSE)

##################NPK_root
orders <- unique(Barley_df_filtered_root_NPK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Barley_df_filtered_root_NPK$Order)) {
  idx <- Barley_df_filtered_root_NPK$Order %in% c(order)
  Barley_df_filtered <- Barley_df_filtered_root_NPK[idx, ]
  if (nrow(Barley_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Barley_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "BarleyCSSP_root_NPK_ANOVA_comb.xlsx", row.names = FALSE)

##################PK_root
orders <- unique(Barley_df_filtered_root_PK$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")
final_df <- data.frame(Genotype = genotype_order)

for (order in unique(Barley_df_filtered_root_PK$Order)) {
  idx <- Barley_df_filtered_root_PK$Order %in% c(order)
  Barley_df_filtered <- Barley_df_filtered_root_PK[idx, ]
  if (nrow(Barley_df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = Barley_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "BarleyCSSP_root_PK_ANOVA_comb.xlsx", row.names = FALSE)

##################UF_root
orders <- unique(Barley_df_filtered_root_UF$Order)

# Initialize an empty list to store results
anova_results <- list()

generate_label_df <- function(pairwise, variable) {
  if (is.null(pairwise[[variable]])) {
    return(data.frame(Genotype = genotype_order, Letters = NA))  # Return NAs if no results
  }
  
  Tukey.levels <- pairwise[[variable]][, 4]  # Extract p-values
  if (any(is.na(Tukey.levels))) {
    return(data.frame(Genotype = genotype_order, Letters = NA))  # Handle missing values
  }
  
  Tukey.labels <- multcompLetters(Tukey.levels)$Letters  # Get letter groupings
  Tukey.labels.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)  # Convert to df
  return(Tukey.labels.df)
}

# Define genotype order
genotype_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")

# Initialize final dataframe
final_df <- data.frame(Genotype = genotype_order)

for (order in orders) {
  idx <- Barley_df_filtered_root_UF$Order %in% c(order)
  Barley_df_filtered <- Barley_df_filtered_root_UF[idx, ]
  if (nrow(Barley_df_filtered) < length(genotype_order)) {
    next
  }
  ano <- aov(RA ~ Genotype, data = Barley_df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "BarleyCSSP_root_UF_ANOVA_comb.xlsx", row.names = FALSE)

