# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the required packages.
pkg <- c("ggplot2", "dplyr", "forcats", "multcompView")
for(pk in pkg){
  library(pk, character.only = T)
}

# Load the chao1 and metadata file.
alpha <- read.table("../../../2_rarefication_chao1/2_Hordeum/2_chao1/HordeumCSSP_AskovSoils_chao1.txt",
                    sep = "\t",       
                    header = TRUE,    
                    row.names = 1,    
                    check.names = FALSE)

design <- read.table("../../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
                     sep = "\t",       
                     header = TRUE,    
                     row.names = 1,    
                     check.names = FALSE)

# Combine the alpha diversity (chao1) and metadata info in new "index" file.
index <- cbind(alpha[, 1],
               design[match(row.names(alpha), row.names(design)), ])
colnames(index)[1] <- "value"

# Filter to only keep WT data.
index <- index %>% filter(Genotype == "WT")

# Set the soil and compartment factor levels.
l1 <- c("NPK","PK","UF")
index$Soil <- factor(index$Soil, levels=l1)
index$Compartment <- factor(index$Compartment,
                            levels = c("Rhizosphere","Root"))

# Define the colours for the plot.
colors <- data.frame(group=l1, color=c("#6F944F","#B2563C","#3C7D82"))
colors <- colors[match(l1, colors$group), ]

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

## Use the function for each compartment.
for(comp in levels(index$Compartment)) {
  
  sub_df <- subset(index, Compartment == comp)
  
  ## Make a summary dataframe for all outputs.
  summary_df <- sub_df %>%
    group_by(Soil) %>%
    summarise(Mean=mean(value),
              Max=max(value),
              Min=min(value),
              Median=median(value),
              Std=sd(value),
              .groups="drop")
  
  ## Perform the ANOVA and TukeyHSD.
  ano <- aov(value ~ Soil, data=sub_df)
  pairwise <- TukeyHSD(ano)
  letters <- generate_label_df(pairwise, "Soil")
  
  ## Add the letters and define the y-positions.
  # y_offset <- 0.075 * (max(summary_df$Max) - min(summary_df$Min))
  y_offset <- 125
  
  tmp <- summary_df %>%
    mutate(Letters = letters$Letters,
           y_position = Max + y_offset,
           Compartment = comp)
  
  label_df <- rbind(label_df, tmp)
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

# Make the plot.
p1 <- ggplot(index, aes(x=Soil, y=value, fill=Soil)) +
  geom_boxplot(alpha=0.7, position=position_dodge(width = 0.7), outlier.color=NA, width=0.3) +
  geom_jitter(position=position_jitter(width=0, height=0.17), size=1, alpha=1) +
  scale_fill_manual(values=as.character(colors$color)) +
  labs(x="", y="Chao1 index", title="Hordeum") +
  geom_text(data=label_df, aes(x=Soil, y=y_position, label=Letters), inherit.aes=FALSE, size=8/.pt) +
  facet_wrap(vars(fct_relevel(Compartment, "Rhizosphere", "Root")),
             scales = "fixed") +
  main_theme+
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0))

p1

# Save the plot.
ggsave("Hordeum_AskovSoils_WT_chao1_rfd.pdf", p1, width=7, height=5, units="cm")
saveRDS(p1, file = "Hordeum_AskovSoils_WT_chao1_rfd.rds")
saveRDS(p1, file = "../../7_final_figures/Hordeum_AskovSoils_WT_chao1_rfd.rds")

# Save ANOVA Tukey HSD output file.
write.csv(label_df, file =  "Hordeum_AskovSoils_WT_chao1_ANOVA_TukeyHSD_rfd.csv")

