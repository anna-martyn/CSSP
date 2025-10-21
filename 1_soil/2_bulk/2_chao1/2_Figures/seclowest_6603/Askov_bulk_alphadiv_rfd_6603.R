# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load chao1 and metadata file.
alpha <- read.table("alpha-diversity.tsv",
                    sep = "\t",       
                    header = TRUE,    
                    row.names = 1,    
                    check.names = FALSE)

design <- read.table("BarleyCSSP_Askov_reseq_metadata.txt",
                     sep = "\t",       
                     header = TRUE,    
                     row.names = 1,    
                     check.names = FALSE)

# Combine alpha diversity (chao1) info and metadata info in new "index" file.
index <- cbind(alpha[, 1],
               design[match(row.names(alpha), row.names(design)), ])
colnames(index)[1] <- "value"

# Filter to only keep soil data and remove Lotus bulk data.
index_bulk <- subset(index,
                     Genotype == "Soil" &
                       !grepl("Lj", Description, ignore.case = FALSE))

# Load package for plotting.
library(ggplot2)

# Define colours and groups for plot.
# colors <- data.frame(group=c("NPK","PK","UF"), color=c("#6D3B00","#A06A37","#D2B48C"))
# colors <- data.frame(group=c("NPK","PK","UF"), color=c("#341C02","#A06A37","#D2B48C"))
# colors <- data.frame(group=c("NPK","PK","UF"),
#                      color=c("#73675A","#A89F90","#D6C9B8"))
# colors <- data.frame(group=c("NPK","PK","UF"),
#                      color=c("#9CAF88","#D08970","#7A9E9F"))
# colors <- data.frame(group=c("NPK","PK","UF"),
#                      colors=c("#A8D17F","#E07B5F","#5FB0B7"))
colors <- data.frame(group=c("NPK","PK","UF"),
                     colors=c("#6F944F","#B2563C","#3C7D82"))
l1 <- c("NPK", "PK", "UF")
index_bulk$Soil <- factor(index_bulk$Soil, levels=l1)
colors <- colors[match(l1, colors$group), ]

# Make plot.
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
                    legend.background=element_blank(),
                    plot.title = element_text(size=8, hjust=1))

# Load package for statistics.
library(dplyr)
library(car)
library(multcompView)

# Make summary file of minimum and maximum values, mean, median and standard deviation.
chao_summarize <- index_bulk%>%
  group_by(Soil)%>%
  summarise(Mean=mean(value), Max=max(value), Min=min(value), Median=median(value), Std=sd(value))

## ANOVA and Tukey HSD.
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
label_df <- chao_summarize %>%
  mutate(Letters = labels$Letters,
         y_position = Max + 0.075 * (max(Max) - min(Min))) %>%
  select(Soil, Max, Letters, y_position)

# Extract overall ANOVA p-value and prepare for plotting.
anova_p <- summary(ano)[[1]][["Pr(>F)"]][1]

# Make plot with significance letters and ANOVA p-value as plot title.
p1 <- ggplot(index_bulk, aes(x=Soil, y=value, fill=Soil)) +
  geom_boxplot(alpha=0.7, position=position_dodge(width = 0.7),
               outlier.color=NA, width=0.3) +
  geom_jitter(position=position_jitter(width=0, height=0.17),
              size=1, alpha=1) +
  scale_fill_manual(values=as.character(colors$color)) +
  labs(x="", y="Chao1 index",
       title = paste0("p = ", signif(anova_p, 3))) +  
  # geom_text(data=label_df, aes(x=Soil, y=y_position, label=Letters),
  #           inherit.aes=FALSE, size=6) +
  main_theme


p1

# Save the plot.
ggsave("Askov_barley_bulk_chao1_rfd_6603.pdf", 
       p1, width = 3, height = 5, units = "cm")
saveRDS(p1, file = "Askov_barley_bulk_chao1_rfd_6603.rds")
saveRDS(p1, file = "../../../6_final_figure/Askov_barley_bulk_chao1_rfd_6603.rds")

# Save ANOVA Tukey HSD output file.
write.csv(labels, file =  "Bulk_chao1_ANOVA_TukeyHSD_rfd_6603.csv")

