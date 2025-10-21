# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load files.
design.file <- paste("Lotus_CSSP_AskovSoils_metadata_excl_new_bulkUF.txt", sep="")
taxonomy.file <- paste("LotusSep_exclUFnew_10_4_silva138_taxonomy-3.txt", sep="")
otu_table.file <- paste("LotusSep_exclUFnew_ASVtable_10_4_nospike.txt", sep="")

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

# Load required packages.
library(dplyr)
library(data.table)
library(ggplot2)
library(multcompView)
library(openxlsx)

# Clean-up taxonomy file layout.
taxa <- lapply(as.list(taxonomy$Taxon),
               function(x) x %>% strsplit(split = "; ") %>% unlist() %>% substr(start = 4, stop = 100))

taxa <- lapply(taxa, function(x) c( x, rep("Unknown", 7 - length(x)) ))
taxa <- lapply(taxa, function(x) 
  data.table(Kingdom = x[1],
             Phylum = x[2],
             Class = x[3],
             Order = x[4],
             Family = x[5],
             Genus = x[6],
             Species = x[7])
)
taxa <- rbindlist(taxa)
taxonomy <- data.table(ASV = taxonomy$`Feature.ID`, taxa, Confidence = taxonomy$Confidence)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- rownames(otu_table) %in% taxonomy$ASV
otu_table <- otu_table[idx, ]

idx <- match(rownames(otu_table), taxonomy$ASV)
taxonomy <- taxonomy[idx, ]

# otu_table normalization
otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))
df <- as.data.frame(otu_table_norm)
df$ASV <- row.names(df)
row.names(df) <- NULL

### reshape data: the aim of this step is to transfer wide data to long data(check online the difference of wide and long dataset in r)
df.long <- reshape(df, varying = 1:204, direction = "long", idvar = 'ASV', timevar="sampleID", v.names=c("RA"), 
                   times=c(colnames(df[,1:204])), sep="")
row.names(df.long) <- NULL
df.long <- as.data.frame(df.long)

# Add taxonomy information (order).
df.long_order <- merge(df.long, taxonomy[, c("ASV", "Order")], by = "ASV", all.x = TRUE)

###  make a dataframe that sum RA by tax and sampleID
df.long_order <- df.long_order %>% select (2:4) ### select the column 2,3,4 from df.long dataframe

df.long_order <- df.long_order %>%
  group_by_(.dots = c("Order", "sampleID")) %>%   ### this is to summarize RA by both tax and sampleID
  summarise_all(funs(sum))

### Add genotype information for df.long.
df3 <- data.frame(sampleID = design$SampleID, Genotype=design$Genotype, Compartment = design$Compartment, Plant=design$Plant, Soil=design$Soil)

df.long_order <- full_join(df.long_order, df3, by="sampleID")


idx <- df.long_order$Genotype%in% c("WT","symrk","ccamk","nsp1","nsp2")
df.long_order <- df.long_order[idx,]

idx <- df.long_order$Soil%in% c("NPK")
df.long_order_NPK <- df.long_order[idx,]

idx <- df.long_order_NPK$Compartment%in% c("rhizosphere")
df.long_order_NPK_rhizo <- df.long_order_NPK[idx,]

idx <- df.long_order_NPK$Compartment%in% c("root")
df.long_order_NPK_root <- df.long_order_NPK[idx,]


idx <- df.long_order$Soil%in% c("PK")
df.long_order_PK <- df.long_order[idx,]

idx <- df.long_order_PK$Compartment%in% c("rhizosphere")
df.long_order_PK_rhizo <- df.long_order_PK[idx,]

idx <- df.long_order_PK$Compartment%in% c("root")
df.long_order_PK_root <- df.long_order_PK[idx,]

idx <- df.long_order$Soil%in% c("UF")
df.long_order_UF <- df.long_order[idx,]

idx <- df.long_order_UF$Compartment%in% c("rhizosphere")
df.long_order_UF_rhizo <- df.long_order_UF[idx,]

idx <- df.long_order_UF$Compartment%in% c("root")
df.long_order_UF_root <- df.long_order_UF[idx,]

# Caluclate mean relative abundance of all orders taking all genotypes into account in the individual soil-fraction combinations.
# Then check which ones have a mean relative abundance >1% in any of the conditions.

## Calculate the mean RA.
process_df <- function(df) {
  df %>%
    group_by(Order) %>%
    summarize(Mean = mean(RA, na.rm = TRUE)) %>%
    arrange(desc(Mean)) %>%
    filter(Mean >= 0.01) %>%
    as.data.frame()
}

# Apply the function to each dataframe
mean_RA_order_NPK_rhizo_0.01 <- process_df(df.long_order_NPK_rhizo)
mean_RA_order_PK_rhizo_0.01  <- process_df(df.long_order_PK_rhizo)
mean_RA_order_UF_rhizo_0.01  <- process_df(df.long_order_UF_rhizo)
mean_RA_order_NPK_root_0.01  <- process_df(df.long_order_NPK_root)
mean_RA_order_PK_root_0.01   <- process_df(df.long_order_PK_root)
mean_RA_order_UF_root_0.01   <- process_df(df.long_order_UF_root)

## Check all unique orders above 1% mean relative abundance in any of the created dataframes.
unique_orders <- bind_rows(mean_RA_order_NPK_rhizo_0.01, mean_RA_order_PK_rhizo_0.01, mean_RA_order_UF_rhizo_0.01, mean_RA_order_NPK_root_0.01, mean_RA_order_PK_root_0.01, mean_RA_order_UF_root_0.01) %>%
  distinct(Order) %>%
  pull(Order)

# Now filter the original dataframe to only keep the orders of interest (RA >1%),a nd only focus on rhizosphere and root fractions.
df_filtered <- df.long_order[df.long_order$Order %in% unique_orders, ]

df_filtered <- df_filtered[df_filtered$Compartment %in% c("rhizosphere", "root"), ]

df_filtered$Genotype <- factor(df_filtered$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))
# df_filtered$Order <- factor(df_filtered$Order, levels = sort(unique(df_filtered$Order)))
df_filtered$Order <- factor(df_filtered$Order, levels = rev(sort(unique(df_filtered$Order))))

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

p1 <- ggplot(df_filtered, aes(x = Genotype, y = Order, fill = RA)) +
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

ggsave(paste("LotusCSSP_heatmap_orders_above_equal1perc.png", sep=""), p1, width=14, height=12)
ggsave(paste("LotusCSSP_heatmap_orders_above_equal1perc.pdf", sep=""), p1, width=14, height=12)

##########SIGNIFICANCE TESTING##################
df_filtered_rhizo <- df_filtered[df_filtered$Compartment %in% c("rhizosphere"), ]
df_filtered_root <- df_filtered[df_filtered$Compartment %in% c("root"), ]

df_filtered_rhizo_NPK <- df_filtered_rhizo[df_filtered_rhizo$Soil %in% c("NPK"), ]
df_filtered_rhizo_PK <- df_filtered_rhizo[df_filtered_rhizo$Soil %in% c("PK"), ]
df_filtered_rhizo_UF <- df_filtered_rhizo[df_filtered_rhizo$Soil %in% c("UF"), ]

df_filtered_root_NPK <- df_filtered_root[df_filtered_root$Soil %in% c("NPK"), ]
df_filtered_root_PK <- df_filtered_root[df_filtered_root$Soil %in% c("PK"), ]
df_filtered_root_UF <- df_filtered_root[df_filtered_root$Soil %in% c("UF"), ]


##################NPK_rhizo
orders <- unique(df_filtered_rhizo_NPK$Order)

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

for (order in unique(df_filtered_rhizo_NPK$Order)) {
  idx <- df_filtered_rhizo_NPK$Order %in% c(order)
  df_filtered <- df_filtered_rhizo_NPK[idx, ]
  if (nrow(df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_rhizo_NPK_ANOVA.xlsx", row.names = FALSE)

##################PK_rhizo
orders <- unique(df_filtered_rhizo_PK$Order)

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

for (order in unique(df_filtered_rhizo_PK$Order)) {
  idx <- df_filtered_rhizo_PK$Order %in% c(order)
  df_filtered <- df_filtered_rhizo_PK[idx, ]
  if (nrow(df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_rhizo_PK_ANOVA.xlsx", row.names = FALSE)

##################UF_rhizo
orders <- unique(df_filtered_rhizo_UF$Order)

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

for (order in unique(df_filtered_rhizo_UF$Order)) {
  idx <- df_filtered_rhizo_UF$Order %in% c(order)
  df_filtered <- df_filtered_rhizo_UF[idx, ]
  if (nrow(df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_rhizo_UF_ANOVA.xlsx", row.names = FALSE)

##################NPK_root
orders <- unique(df_filtered_root_NPK$Order)

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

for (order in unique(df_filtered_root_NPK$Order)) {
  idx <- df_filtered_root_NPK$Order %in% c(order)
  df_filtered <- df_filtered_root_NPK[idx, ]
  if (nrow(df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_root_NPK_ANOVA.xlsx", row.names = FALSE)

##################PK_root
orders <- unique(df_filtered_root_PK$Order)

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

for (order in unique(df_filtered_root_PK$Order)) {
  idx <- df_filtered_root_PK$Order %in% c(order)
  df_filtered <- df_filtered_root_PK[idx, ]
  if (nrow(df_filtered) < length(genotype_order)) {
    next  # Skip if not enough samples
  }
  ano <- aov(RA ~ Genotype, data = df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_root_PK_ANOVA.xlsx", row.names = FALSE)

##################UF_root
orders <- unique(df_filtered_root_UF$Order)

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
  idx <- df_filtered_root_UF$Order %in% c(order)
  df_filtered <- df_filtered_root_UF[idx, ]
  if (nrow(df_filtered) < length(genotype_order)) {
    next  
  }
  ano <- aov(RA ~ Genotype, data = df_filtered)
  pairwise <- TukeyHSD(ano)
  labels_df <- generate_label_df(pairwise, "Genotype")
  labels_df <- labels_df %>%
    filter(Genotype %in% genotype_order) %>%
    arrange(match(Genotype, genotype_order)) %>%
    select(Letters)
  colnames(labels_df) <- order
  final_df <- cbind(final_df, labels_df)
}
write.xlsx(final_df, file = "LotusCSSP_root_UF_ANOVA.xlsx", row.names = FALSE)
