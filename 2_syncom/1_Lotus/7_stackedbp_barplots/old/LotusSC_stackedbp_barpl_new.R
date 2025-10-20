options(warn=-1)

# cleanup
rm(list=ls())

# load packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggh4x)
library(stringr)
library(multcompView)

# directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# laod files
design.file <- paste("LotusCSSP_LjSC_metadata.txt", sep="")
taxonomy.file <- paste("LjSC_taxonomy.txt", sep="")
otu_table.file <- paste("LotusCSSP_LjSC_ASVtable.txt", sep="")

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, row.names=1, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

idx <- rownames(otu_table) %in% taxonomy[,1]
otu_table <- otu_table[idx, ]

idx <- match(rownames(otu_table), taxonomy[,1])
taxonomy <- taxonomy[idx, ]

# otu_table normalization
otu_table_norm <- apply(otu_table, 2, function(x) x / sum(x))
df <- as.data.frame(otu_table_norm)

# add taxonomy info to "df"
tax_order <- taxonomy[ ,5]
df$order <- tax_order

ASVid <- row.names(df)
df <- cbind(ASVid, df)
row.names(df) <- NULL

# reshape data: the aim of this step is to transfer wide data to long data(check online the difference of wide and long dataset in r)
df.long <- reshape(df, varying = 2:60, direction = "long", idvar = 'ASVid', timevar="sampleID", v.names=c("RA"), 
times=c(colnames(df[,2:60])), sep="")
row.names(df.long) <- NULL
df.long <- as.data.frame(df.long)

# make a dataframe that sums RA by taxonomic order and sampleID
df.long_order <- df.long %>% select (2:4) ### select the column 2,3,4 from df.long dataframe

df.long_order <- df.long_order %>%
  group_by_(.dots = c("order", "sampleID")) %>%   ### this is to summarize RA by both tax and sampleID
  summarise_all(funs(sum))

# Add genotype information to df.long
df2 <- data.frame(sampleID = design$SampleID, Genotype=design$Genotype, Compartment = design$Compartment)
df.long_order <- full_join(df.long_order, df2, by="sampleID")

# Next we pick the genotypes we want to focus on, as well as make separate dataframes for rhizosphere and root data.
idx <- df.long_order$Genotype%in% c("WT","symrk","ccamk","nsp1","nsp2")
df.long_order_nof6h1 <- df.long_order[idx,]

idx <- df.long_order_nof6h1$Compartment%in% c("rhizosphere")
df.long_order_rhizo <- df.long_order_nof6h1[idx,]

idx <- df.long_order_nof6h1$Compartment%in% c("root")
df.long_order_endo <- df.long_order_nof6h1[idx,]

idx <- df.long_order_nof6h1$Compartment%in% c("nodules")
df.long_order_nod <- df.long_order_nof6h1[idx,]

## select the top 10 taxa by relative abundance for the stack bar plot

# calculate the mean RA for each order across samples
mean_RA_order_nof6h1 <- df.long_order_nof6h1 %>%
  group_by(order) %>%
  dplyr::summarize(Mean = mean(RA, na.rm = TRUE))

mean_RA_order_nof6h1 <- as.data.frame(mean_RA_order_nof6h1)

mean_RA_order_nof6h1 <- mean_RA_order_nof6h1[order(-mean_RA_order_nof6h1$Mean),]

# Check top 20 orders with highest abundance
mean_RA_order_nof6h1 %>% top_n(20)

# Check how many unique orders in general
unique(mean_RA_order_nof6h1$order)

# stacked bar plot for the top orders (here we simply choose all 9 orders present)

## Give genotype and compartment order for plots.
df.long_order_nof6h1$Genotype <- factor(df.long_order_nof6h1$Genotype, levels = c("WT","symrk","ccamk","nsp1","nsp2"))
df.long_order_nof6h1$Compartment <- factor(df.long_order_nof6h1$Compartment, levels = c("rhizosphere","root","nodules"))

## Define colours for orders.
colors <- data.frame(group=c("Actinomycetales","Bacillales","Burkholderiales","Caulobacterales","Flavobacteriales","Pseudomonadales","Rhizobiales","Sphingomonadales","Xanthomonadales"), 
                      colors=c("#2A0134","#771155", "#645394", "#AA4488","#ffeeef", "#88CCAA", "lightyellow", "#fdbb6b", "#ffc0cb")) 
colors <- colors[colors$group %in%df.long_order_nof6h1$order, ]

## make plot
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans"))
p1 <- ggplot(df.long_order_nof6h1, aes(x=sampleID, y = RA, fill = order)) +
  geom_bar(stat = "identity", width = .5) +
  facet_nested(~Compartment + Genotype,scales ="free_x", labeller = label_parsed)+
  scale_fill_manual(values=colors$colors)+
  main_theme+
  ylab("Relative abundance")+
  theme(legend.position = "bottom")+ guides(fill=guide_legend())+ 
  theme(axis.text.x = element_blank(),
        plot.title = element_text(size = 20, face = "bold"), 
        strip.text.x = element_text(size = 20, face = "bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20))
p1

ggsave(paste("LotusSC_stacked_barplot.png", sep=""), p1, width=14, height=8)
ggsave(paste("LotusSC_stacked_barplot.pdf", sep=""), p1, width=14, height=8)

# boxplots orders

## filter nodule samples out
idx <- df.long_order_nof6h1$Compartment%in% c("rhizosphere","root")
df.long_order_nonod <- df.long_order_nof6h1[idx,]

## Define colours for the genotypes.
colors2 <- data.frame(group=c("WT", "symrk","ccamk","nsp1", "nsp2"), 
                     color=c("#A9C289","#FEDA8B","#FDB366","#C0E4EF","#6EA6CD"))

# colors2 <- data.frame(group=c("WT","symrk","ccamk","nsp1","nsp2"), 
#                      color=c("#33a02c","#ff7f00","#1f78b4","#e31a1c", "#ffd700"))

## make the plot
dodge <- position_dodge (width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.text = element_text(size=20, color = "black"),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans",size=20, color="black"),
                    axis.title.y = element_text(size = 20))

p2 <- ggplot(df.long_order_nonod, aes(x=order, y=RA, fill=Genotype)) +
  geom_boxplot( width=0.5,position = dodge, outlier.color = NA)+
  geom_jitter(aes(group = Genotype), position = position_dodge(width = 0.9), size = 1, alpha = 0.5)+
  scale_fill_manual(values=as.character(colors2$color)) +
  labs(x="", y="Relative Abundance") +
  scale_y_continuous(limits = c(0, 0.7))+
  facet_wrap(~Compartment,  ncol=1)+
  main_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20, angle = 50, hjust = 1))
p2

ggsave(paste("LotusSC_boxplots.png", sep=""), p2, width=16, height=12)
ggsave(paste("LotusSC_boxplots.pdf", sep=""), p2, width=16, height=12)

# barplots

##rhizosphere
idx <- df.long_order_nof6h1$Compartment%in% c("rhizosphere")
df.long_order_rhizo <- df.long_order_nof6h1[idx,]

### Calculate standard deviation and standard error for plots.
df.long_order_summary <- df.long_order_rhizo %>%
  group_by(order, Genotype) %>%
  summarise(
    Mean_RA = mean(RA, na.rm = TRUE),
    SD_RA = sd(RA, na.rm = TRUE),
    SE_RA = SD_RA / sqrt(n())
  )

### plot
dodge <- position_dodge (width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.text = element_text(size=20, color = "black"),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans",size=20, color="black"),
                    axis.title.y = element_text(size = 20))

p3 <- ggplot() +
  geom_bar(data = df.long_order_summary, aes(x=order, y=Mean_RA, fill=Genotype),
           stat = "identity", position = position_dodge(width = 0.9), width = 0.8, alpha=0.8) +
  geom_errorbar(data = df.long_order_summary, 
                aes(x=order, ymin = Mean_RA - SE_RA, ymax = Mean_RA + SE_RA, group=Genotype),
                width = 0.3, position = position_dodge(width = 0.9), color = "black") +
  geom_point(data = df.long_order_rhizo, 
             aes(x=order, y=RA, fill=Genotype),
             position = position_dodge(width = 0.9),
             size =2, alpha = 0.5, shape = 21, color = "black") +  
  scale_fill_manual(values=as.character(colors2$color)) +
  scale_color_manual(values=as.character(colors2$color)) +
  labs(x="", y="Relative Abundance") +
  scale_y_continuous(limits = c(0, 0.6)) +
  main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 70, hjust=1))

p3

ggsave(paste("LotusSC_boxplots_rhizo.png", sep=""), p3, width=20, height=10)
ggsave(paste("LotusSC_boxplots_rhizo.pdf", sep=""), p3, width=20, height=10)

##root
idx <- df.long_order_nof6h1$Compartment%in% c("root")
df.long_order_root <- df.long_order_nof6h1[idx,]

### Calculate standard deviation and standard error for plots.
df.long_order_summary <- df.long_order_root %>%
  group_by(order, Genotype) %>%
  summarise(
    Mean_RA = mean(RA, na.rm = TRUE),
    SD_RA = sd(RA, na.rm = TRUE),
    SE_RA = SD_RA / sqrt(n())
  )

### plot
dodge <- position_dodge (width = 0.9)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.grid=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=20),
                    legend.text = element_text(size=20, color = "black"),
                    legend.position="right",
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans",size=20, color="black"),
                    axis.title.y = element_text(size = 20))

p3 <- ggplot() +
  geom_bar(data = df.long_order_summary, aes(x=order, y=Mean_RA, fill=Genotype),
           stat = "identity", position = position_dodge(width = 0.9), width = 0.8, alpha=0.8) +
  geom_errorbar(data = df.long_order_summary, 
                aes(x=order, ymin = Mean_RA - SE_RA, ymax = Mean_RA + SE_RA, group=Genotype),
                width = 0.3, position = position_dodge(width = 0.9), color = "black") +
  geom_point(data = df.long_order_root, 
             aes(x=order, y=RA, fill=Genotype),
             position = position_dodge(width = 0.9),
             size =2, alpha = 0.5, shape = 21, color = "black") +  
  scale_fill_manual(values=as.character(colors2$color)) +
  scale_color_manual(values=as.character(colors2$color)) +
  labs(x="", y="Relative Abundance") +
  scale_y_continuous(limits = c(0, 0.7)) +
  main_theme +
  theme(legend.position= "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 70, hjust=1))

p3

ggsave(paste("LotusSC_boxplots_root.png", sep=""), p3, width=20, height=10)
ggsave(paste("LotusSC_boxplots_root.pdf", sep=""), p3, width=20, height=10)

## Statistical analysis

###Rhizosphere
# Get unique orders from the dataframe
unique_orders <- unique(df.long_order_rhizo$order)

# Initialize an empty list to store results
results_list <- list()

# Loop through each order and perform ANOVA
for (order_name in unique_orders) {
  
  # Subset dataframe for the current order
  df_subset <- df.long_order_rhizo[df.long_order_rhizo$order == order_name, ]
  
  # Check if there are multiple Genotypes to compare
  if (length(unique(df_subset$Genotype)) > 1) {
    
    # Perform ANOVA
    ano <- aov(RA ~ Genotype, data = df_subset)
    
    # Check if ANOVA is significant before proceeding
    if (anova(ano)[["Pr(>F)"]][1] < 0.05) {
      
      # Perform Tukey's HSD test
      pairwise <- TukeyHSD(ano)
      
      # Function to generate a single letter per genotype group
      generate_label_df <- function(pairwise, variable) {
        Tukey.levels <- pairwise[[variable]][, 4]  # P-values
        if (length(Tukey.levels) > 0) {
          Tukey.labels <- multcompLetters(Tukey.levels)$Letters
          Tukey.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)
          return(Tukey.df)
        } else {
          return(NULL)
        }
      }
      
      LABELS <- generate_label_df(pairwise, "Genotype")
      
      if (!is.null(LABELS)) {
        colnames(LABELS)[2] <- order_name  # Rename column to order name
        results_list[[order_name]] <- LABELS
      } else {
        cat("Skipping", order_name, "- No valid Tukey test results.\n")
      }
      
    } else {
      cat("Skipping", order_name, "- ANOVA not significant.\n")
    }
    
  } else {
    cat("Skipping", order_name, "- Not enough genotypes.\n")
  }
}

# Define the desired genotype order
desired_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")

# Merge results only if there is valid data
if (length(results_list) > 0) {
  final_results <- Reduce(function(x, y) merge(x, y, by = "Genotype", all = TRUE), results_list)
  
  # Reorder genotypes based on the predefined order
  final_results <- final_results[match(desired_order, final_results$Genotype), ]
  
  # Save to CSV
  write.csv(final_results, file = "LotusSC_rhizo_ANOVA_summary.csv", row.names = FALSE)
  
  print("ANOVA summary file saved: LotusSC_rhizo_ANOVA_summary.csv")
} else {
  print("No significant ANOVA results found. No file was saved.")
}


###Root
# Get unique orders from the dataframe
unique_orders <- unique(df.long_order_root$order)

# Initialize an empty list to store results
results_list <- list()

# Loop through each order and perform ANOVA
for (order_name in unique_orders) {
  
  # Subset dataframe for the current order
  df_subset <- df.long_order_root[df.long_order_root$order == order_name, ]
  
  # Check if there are multiple Genotypes to compare
  if (length(unique(df_subset$Genotype)) > 1) {
    
    # Perform ANOVA
    ano <- aov(RA ~ Genotype, data = df_subset)
    
    # Extract p-value and check if it's NA
    p_value <- anova(ano)[["Pr(>F)"]][1]
    
    if (!is.na(p_value) && p_value < 0.05) {
      
      # Perform Tukey's HSD test
      pairwise <- TukeyHSD(ano)
      
      # Function to generate a single letter per genotype group
      generate_label_df <- function(pairwise, variable) {
        Tukey.levels <- pairwise[[variable]][, 4]  # P-values
        if (length(Tukey.levels) > 0) {
          Tukey.labels <- multcompLetters(Tukey.levels)$Letters
          Tukey.df <- data.frame(Genotype = names(Tukey.labels), Letters = Tukey.labels)
          return(Tukey.df)
        } else {
          return(NULL)
        }
      }
      
      LABELS <- generate_label_df(pairwise, "Genotype")
      
      if (!is.null(LABELS)) {
        colnames(LABELS)[2] <- order_name  # Rename column to order name
        results_list[[order_name]] <- LABELS
      } else {
        cat("Skipping", order_name, "- No valid Tukey test results.\n")
      }
      
    } else {
      cat("Skipping", order_name, "- ANOVA not significant or p-value is NA.\n")
    }
    
  } else {
    cat("Skipping", order_name, "- Not enough genotypes.\n")
  }
}

# Define the desired genotype order
desired_order <- c("WT", "symrk", "ccamk", "nsp1", "nsp2")

# Merge results only if there is valid data
if (length(results_list) > 0) {
  final_results <- Reduce(function(x, y) merge(x, y, by = "Genotype", all = TRUE), results_list)
  
  # Reorder genotypes based on the predefined order
  final_results <- final_results[match(desired_order, final_results$Genotype), ]
  
  # Save to CSV
  write.csv(final_results, file = "LotusSC_root_ANOVA_summary.csv", row.names = FALSE)
  
  print("ANOVA summary file saved: LotusSC_root_ANOVA_summary.csv")
} else {
  print("No significant ANOVA results found. No file was saved.")
}