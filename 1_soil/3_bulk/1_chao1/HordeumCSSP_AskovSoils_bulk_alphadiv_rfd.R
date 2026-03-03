# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory to source file location.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load required packages.
library(ggplot2)
library(dplyr)
library(car)
library(multcompView)

# Load chao1 and metadata file.
alpha <- read.table("../../2_rarefication_chao1/2_Hordeum/2_chao1/HordeumCSSP_AskovSoils_chao1.txt",
                    sep = "\t",       
                    header = TRUE,    
                    row.names = 1,    
                    check.names = FALSE)

design <- read.table("../../1_data/2_Hordeum/HordeumCSSP_AskovSoils_metadata.txt",
                     sep = "\t",       
                     header = TRUE,    
                     row.names = 1,    
                     check.names = FALSE)  

# Combine alpha diversity (chao1) info and metadata info in new "index" file.
index <- cbind(alpha[, 1],
               design[match(row.names(alpha), row.names(design)), ])
colnames(index)[1] <- "value"

# Filter to only keep soil data.
index_bulk <- subset(index,
                     Genotype == "Soil")

# Define colours and groups for plot.
colors <- data.frame(group=c("NPK","PK","UF"), color=c("#6F944F","#B2563C","#3C7D82"))
l1 <- c("NPK", "PK", "UF")
index_bulk$Soil <- factor(index_bulk$Soil, levels=l1)
colors <- colors[match(l1, colors$group), ]

# Set main theme for the plot.
main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 6, color = "black"),
                    legend.text = element_text(size = 6, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 6),
                    legend.position="none",
                    legend.background=element_blank(),
                    plot.title = element_text(size = 6, hjust=1))

# Make summary file of minimum and maximum values, mean, median and standard deviation.
chao_summarise <- index_bulk%>%
  group_by(Soil)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))

# Perfom significance analysis using ANOVA and Tukey HSD.
ano <- aov(value ~ Soil, data=index_bulk)
anova(ano)
pairwise <- TukeyHSD(ano)
generate_label_df <- function(pairwise, variable){
  Tukey.levels <- pairwise[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey.labels$Type=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$Type) , ]
  return(Tukey.labels)
}
labels=generate_label_df(pairwise , "Soil")

# Prepare label positions above the max value per group.
label_df <- chao_summarise %>%
  mutate(Letters = labels$Letters,
         y_position = Max + 0.075 * (max(Max) - min(Min))) %>%
  select(Soil, Max, Letters, y_position)

# Extract overall ANOVA p-value and prepare for plotting.
anova_p <- summary(ano)[[1]][["Pr(>F)"]][1]

# Make plot with significance letters and ANOVA p-value as plot title.
p1 <- ggplot(index_bulk, aes(x=Soil, y=value, fill=Soil)) +
  geom_boxplot(alpha=0.7, position=position_dodge(width = 0.7), outlier.color=NA, width=0.3) +
  geom_jitter(position=position_jitter(width=0, height=0.17), size=1, alpha=1) +
  scale_fill_manual(values=as.character(colors$color)) +
  labs(x="", y="Chao1 index") + 
  geom_text(data=label_df, aes(x=Soil, y=y_position, label=Letters), inherit.aes=FALSE, size = 6/.pt) +
  main_theme+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  NULL

p1

# Save the plot.
ggsave("HordeumCSSP_AskovSoils_bulk_chao1_rfd.pdf", p1, width=3, height=5, units="cm")
saveRDS(p1, file = "HordeumCSSP_AskovSoils_bulk_chao1_rfd.rds")
saveRDS(p1, file = "../5_final_figure/HordeumCSSP_AskovSoils_bulk_chao1_rfd.rds")

# Save ANOVA and Tukey HSD output file.
write.csv(labels, file =  "HordeumCSSP_AskovSoils_bulk_chao1_ANOVA_TukeyHSD_rfd.csv")

