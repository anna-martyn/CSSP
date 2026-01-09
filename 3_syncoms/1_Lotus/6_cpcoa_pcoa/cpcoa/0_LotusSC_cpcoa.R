# Clean up.
options(warn=-1)
rm(list=ls())

# Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load plotting functions.
library("ggplot2")
library("vegan")
library("ggtext")

# Load cpcoa functions file.
source("cpcoa.func.R")

# Load input files.
design <- read.table("../../1_data/LotusSC_metadata.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)
asv_table <- read.table(
  "../../1_data/LotusSC_ASVtable_nospike.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1
)

# # Keep samples that are present in design and ASV table files for each dataset.
# idx<- intersect(design$SampleID, colnames(asv_table))
# design <- design[match(idx, design$SampleID), ]
# asv_table <- asv_table[, idx]

# Subset datasets to only keep samples for genotypes of interest.
design$SampleID <- row.names(design)
# samples_keep <- design$SampleID[design$Genotype %in% c("WT","symrk","ccamk","nsp1","nsp2")]
# design <- design[design$SampleID %in% samples_keep, ]
# asv_table <- asv_table[, samples_keep]

# Change compartment names in design files.
# design$Compartment <- recode(design$Compartment,
#                              "rhizo"="Rhizosphere",
#                              "endo"="Root",
#                              "nod"="Nodules")

# For plotting the beta-diversity based on Bray Curtis dissimilarities, we will use to approaches:
## 1. Take all ASVs into account (matched ASVs and contaminants).
## 2. Take only matched ASVs into account.

# Make filtered asv_table for approach where only matched ASVs are taken into account.
asv_table_filt <- asv_table[grepl("Lj", row.names(asv_table)), ]

# Set colours and shapes for plot, as well as plotting parameters.
colors <- c("WT"="#A9C289","symrk"="#FEDA8B","ccamk"="#FDB366","nsp1"="#C0E4EF","nsp2"="#6EA6CD")
shapes <- c("Rhizosphere"=15,"Root"=16,"Nodules"=17)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid.major = element_line(color = "gray90"),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                    axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text = element_text(size = 20, color = "black"),
                    legend.text = element_text(size=20, color = "black"),
                    legend.key=element_blank(),
                    axis.title.y = element_text(size = 20),
                    text=element_text(size=20, color="black"),
                    legend.position="none",
                    # legend.background=element_rect(colour="black", fill=NA),
                    legend.background=element_blank(),
                    plot.title = element_text(size=20, hjust=0.9))

# Make cpcoa and plotting functions.
plot_cpcoa <- function(asv_table_norm, design, plot_title){

  # Bray-Curtis distance already implied in capscale
  cap <- capscale(t(asv_table_norm) ~ Compartment*Genotype + Condition(Bio_rep),
                  data=design, add=F, sqrt.dist=T, distance="bray")

  # Permutation ANOVA
  perm_anova <- anova.cca(cap)

  # Variance explained
  var_tbl <- variability_table(cap)
  eig <- cap$CCA$eig
  variance <- var_tbl["constrained","proportion"]
  p.val <- perm_anova[1,4]

  # Weighted average (sample scores)
  points <- as.data.frame(cap$CCA$wa[,1:2])
  colnames(points) <- c("x","y")
  points <- cbind(points, design[match(rownames(points), design$SampleID), ])

  points$Genotype <- factor(points$Genotype, levels=names(colors))
  points$Compartment <- factor(points$Compartment, levels=names(shapes))

  # Centroids
  centroids <- aggregate(cbind(x,y) ~ Genotype + Compartment, data=points, FUN=mean)
  segments <- merge(points, setNames(centroids, c("Genotype","Compartment","seg_x","seg_y")),
                    by=c("Genotype","Compartment"), sort=FALSE)

  # Italic labels for symrk, ccamk, nsp1, nsp2
  genotype_labels <- c(
    "WT"     = "WT",
    "symrk"  = "*symrk*",
    "ccamk"  = "*ccamk*",
    "nsp1"   = "*nsp1*",
    "nsp2"   = "*nsp2*"
  )

  # Plot
  ggplot(points, aes(x=x, y=y, color=Genotype, shape=Compartment)) +
    geom_segment(data=segments, aes(x=x, y=y, xend=seg_x, yend=seg_y, color=Genotype),
                 alpha=0.5, show.legend=FALSE) +
    geom_point(size=5, alpha=0.7) +
    scale_color_manual(values=colors, labels=genotype_labels) +
    scale_shape_manual(values=shapes) +
    labs(x=paste0("CPCo 1 (", format(100*eig[1]/sum(eig), digits=4), "%)"),
         y=paste0("CPCo 2 (", format(100*eig[2]/sum(eig), digits=4), "%)")) +
    ggtitle(plot_title, subtitle=paste0(format(100*variance,digits=3), "% of variance; p=", format(p.val,digits=2))) +
    main_theme +
    guides(color = guide_legend(override.aes = list(size=5)),
           shape = guide_legend(override.aes = list(size=5))) +
    theme(
      axis.title.x = element_text(size=20),
      plot.title = element_text(face="bold", size=20, hjust=0),
      plot.subtitle = element_text(size=20, hjust=0),
      legend.position="right",
      legend.text = element_markdown(size=20)  # <--- ggtext makes mutant names italic
    )
}

# Normalise the ASV tables.
asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))
asv_table_filt_norm <- apply(asv_table_filt, 2, function(x) x/sum(x))

# Execute functions for both datasets.
p_allASVs <- plot_cpcoa(asv_table_norm, design, "All ASVs")
p_filtASVs <- plot_cpcoa(asv_table_filt_norm, design, "Matched ASVs")

# Save plots.
ggsave("LotusSC_cpcoa_allASVs.pdf", p_allASVs, width=5, height=5)
saveRDS(p_allASVs, "LotusSC_cpcoa_allASVs.rds")
saveRDS(p_allASVs, "../../../3_final_figures/LotusSC_cpcoa_allASVs.rds")

ggsave("LotusSC_cpcoa_matchedASVsonly.pdf", p_filtASVs, width=5, height=5)
saveRDS(p_filtASVs, "LotusSC_cpcoa_matchedASVsonly.rds")
saveRDS(p_filtASVs, "../../../3_final_figures/LotusSC_cpcoa_matchedASVsonly.rds")
