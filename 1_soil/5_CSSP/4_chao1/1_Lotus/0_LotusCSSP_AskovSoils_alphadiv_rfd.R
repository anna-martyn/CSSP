# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
library(ggplot2)
library(dplyr)
library(multcompView)
library(forcats)
library(ggtext)

# Load the chao1 and metadata file.
alpha <- read.table("../../../2_rarefication_chao1/1_Lotus/2_chao1/LotusCSSP_AskovSoils_chao1.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
design <- read.table("../../../1_data/1_Lotus/LotusCSSP_AskovSoils_metadata.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)

# Combine the alpha diversity (chao1) and metadata info in new "index" file.
index <- cbind(alpha[,1], design[match(row.names(alpha), row.names(design)), ])
colnames(index)[1] <- "value"

# Set the soil, compartment, and genotype factor levels.
index$Soil <- factor(index$Soil, levels=c("NPK","PK","UF"))
index$Compartment <- factor(index$Compartment, levels=c("Rhizosphere","Root","Nodules"))
index$Genotype <- factor(index$Genotype,
                            levels = c("WT","symrk","ccamk","nsp1","nsp2"))

# Define the colours for the plot.
colors <- data.frame(
  group = c("WT", "symrk","ccamk","nsp1", "nsp2"), 
  color = c("#A9C289", "#FEDA8B", "#FDB366", "#C0E4EF", "#6EA6CD")
)

# Make mutant genotype names italic.
genotype_labels <- c(
  "WT" = "WT",
  "symrk" = "italic(symrk)",
  "ccamk" = "italic(ccamk)",
  "nsp1" = "italic(nsp1)",
  "nsp2" = "italic(nsp2)"
)

# Perform significance analysis using ANOVA and Tukey HSD.

## Define a function for the generation of significance letters.
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$Type = rownames(Tukey.labels)
  Tukey.labels = Tukey.labels[order(Tukey.labels$Type), ]
  return(Tukey.labels)
}

## Prepare an emoty dataframe for the significance labels.
label_df <- data.frame()

## Use the function for each soil-compartment combination.
for(s in levels(index$Soil)) {
  for(c in levels(index$Compartment)) {

    sub_df <- subset(index, Soil == s & Compartment == c)
    
    # Skip if fewer than 2 genotypes present (this is the case for the nodule compartment).
    if(length(unique(sub_df$Genotype)) < 2) {
      next
    }
    
## Make a summary dataframe for all outputs.
  summary_df <- sub_df %>%
    group_by(Genotype) %>%
    summarise(Mean=mean(value),
              Max=max(value),
              Min=min(value),
              Median=median(value),
              Std=sd(value),
              .groups="drop")

  ## Perform the ANOVA and TukeyHSD.
  ano <- aov(value ~ Genotype, data=sub_df)
  pairwise <- TukeyHSD(ano)
  letters <- generate_label_df(pairwise, "Genotype")

  ## Add the letters and define the y-positions.
  y_offset <- 0.3 * (max(summary_df$Max) - min(summary_df$Min))
  
  tmp <- summary_df %>%
    mutate(Letters = letters$Letters,
           y_position = Max + y_offset,
           Soil = s,
           Compartment = c)

  label_df <- rbind(label_df, tmp)
  }
}

# Set the main theme for the plot and make the plot.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 8, color = "black"),
                    legend.text = element_text(size=8, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 8),
                    legend.position="none",
                    strip.text = element_text(size = 8, color = "black"),
                    legend.background=element_blank(),
                    plot.title = element_text(size=8, hjust=1))

p1 <- ggplot(index, aes(x=Genotype, y=value, fill=Genotype)) +
  geom_boxplot(alpha=0.7, position=position_dodge(width = 0.7),
               outlier.color=NA, width=0.3) +
  geom_jitter(position=position_jitter(width=0, height=0.17), size=1, alpha=1) +
  scale_fill_manual(values=as.character(colors$color)) +
  scale_x_discrete(labels = function(x) parse(text = genotype_labels[x])) +
  labs(x="", y="Chao1 index", title="Lotus") +
  geom_text(data=label_df, aes(x=Genotype, y=y_position, label=Letters),
            inherit.aes=FALSE, size=8/.pt) +
  facet_grid(Compartment ~ Soil, scales = "free_y", drop = FALSE) +
  main_theme +
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1))
  

p1

# Save the plot.
ggsave("LotusCSSP_AskovSoils_chao1_rfd.pdf", p1, width=12, height=10, unit="cm")
saveRDS(p1, file = "LotusCSSP_AskovSoils_chao1_rfd.rds")
saveRDS(p1, file = "../../8_final_figures/LotusCSSP_AskovSoils_chao1_rfd.rds")

# Save ANOVA Tukey HSD output file.
write.csv(label_df, file =  "LotusCSSP_AskovSoils_chao1_ANOVA_TukeyHSD_rfd.csv")

