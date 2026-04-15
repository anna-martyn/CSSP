# Loading packages and setting colours -----------------------------------------
pkg <- c(
  "data.table",
  "ggplot2",
  "ggh4x",
  "gridExtra",
  "cowplot",
  "ggtext",
  "htmlTable",
  "Gmisc",
  "grid"
)
for(pk in pkg){
  library(pk, character.only = TRUE)
}

# Setting working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Defining colours for genotypes and orders
colors_genotype <- c(
  "WT" = "#A9C289",
  "symrk" = "#FEDA8B",
  "ccamk" = "#FDB366",
  "nsp1" = "#C0E4EF",
  "nsp2" = "#6EA6CD"
)

colors_order <- fread("../../../0_files/Bacterial_order_colors.csv")

# Loading results from prediction analysis
res <- fread("../1_predictions/1_tables/Prediction_results.csv")
tax_summary <- fread("../1_predictions/1_tables/Pred_taxonomic_composition.csv")
ratio_amount <- fread("../1_predictions/1_tables/Pred_accuracy_summary.csv")

# Flowchart --------------------------------------------------------------------

# Defining a custom function to create a diamond shape, and defining the individual
# boxes and arrows for the flowchart
diamond_grob <- function(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar()){
  diamond_points <- matrix(
    c(0.5, 0, 1, 0.5, 0.5, 1, 0, 0.5),
    ncol = 2,
    byrow = TRUE
  )
  polygonGrob(
    x = diamond_points[, 1],
    y = diamond_points[, 2],
    default.units = "npc",
    gp = gp
  )
}

grid.newpage()

# Making boxes with text
train <- boxGrob(
  "Training\n Data",
  x = 0.55,
  y = 0.9,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

test <- boxGrob(
  label = "Test Data",
  x = 0.28,
  y = 0.9,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

mod1 <- boxGrob(
  label = "Model 1: \nPredicts NPK vs. non-NPK \nTrained on NPK, PK, and UF data",
  x = 0.28,
  y = 0.65,
  width = unit(42, "mm"),
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

mod2 <- boxGrob(
  label = "Model 2: \nPredicts PK vs. UF \nTrained on PK and UF data",
  x = 0.78,
  y = 0.65,
  width = unit(37, "mm"),
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

pred <- boxGrob(
  label = "Prediction on Test Data: \nNPK or non-NPK?",
  x = 0.28,
  y = 0.45,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

# Using diamond_grob as the box_fn in boxGrob
diamond_box <- boxGrob(
  label = "NPK?",
  box_fn = diamond_grob,
  x = 0.28,
  y = 0.3,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

pred2 <- boxGrob(
  label = "Prediction on Test Data: \nPK or UF?",
  x = 0.78,
  y = 0.3,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

concl <- boxGrob(
  label = "Prediction \nconcluded",
  x = 0.78,
  y = 0.1,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

concl2 <- boxGrob(
  label = "Prediction \nconcluded",
  x = 0.28,
  y = 0.1,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "#f5fbff", col = "#4881b3", lwd = 2)
)

no <- boxGrob(
  label = "No",
  x = 0.45,
  y = 0.3,
  width = 0.05,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "white", col = "white")
)

yes <- boxGrob(
  label = "Yes",
  x = 0.28,
  y = 0.22,
  height = 0.025,
  txt_gp = gpar(col = "black", fontsize = 6, family = "helvetica"),
  box_gp = gpar(fill = "white", col = "white")
)

# Printing boxes
train; test; mod1; mod2; pred; diamond_box; pred2; concl; concl2

# Making connections between boxes
## Custom arrow
custom_arrow <- arrow(type = "closed", length = unit(0.3, "cm"))

connect1 <- connectGrob(
  start = train,
  end = mod1,
  "N",
  "l",
  arrow_obj = custom_arrow
)
connect2 <- connectGrob(
  start = train,
  end = mod2,
  "N",
  "l",
  arrow_obj = custom_arrow
)
connect3 <- connectGrob(
  start = mod1,
  end = pred,
  "vertical",
  arrow_obj = custom_arrow
)

connect4 <- connectGrob(start = test, end = pred, "Z", arrow_obj = custom_arrow)
x_coords <- attr(connect4, "line")$x
y_coords <- attr(connect4, "line")$y

# Adjusting end-points of Z-shaped connection
x_coords[2] <- x_coords[2] - unit(1.5, "mm")
x_coords[3] <- x_coords[3] - unit(1.5, "mm")
x_coords[4] <- coords(pred)$left + coords(pred)$width / 2 + unit(1.35, "mm")
y_coords[3] <- y_coords[3] + unit(16, "mm")
y_coords[4] <- y_coords[4] + unit(16, "mm")

new_connection <- linesGrob(
  x = x_coords,
  y = y_coords,
  gp = gpar(fill = "black")
)

connect5 <- connectGrob(
  start = pred,
  end = diamond_box,
  "vertical",
  arrow_obj = custom_arrow
)
connect6 <- connectGrob(
  start = mod2,
  end = pred2,
  "vertical",
  arrow_obj = custom_arrow
)
connect7 <- connectGrob(
  start = diamond_box,
  end = pred2,
  type = "horizontal",
  arrow_obj = custom_arrow,
  lty_gp = gpar(col = "#e81313", fill = "#e81313")
)
connect8 <- connectGrob(
  start = diamond_box,
  end = concl2,
  type = "vertical",
  arrow_obj = custom_arrow,
  lty_gp = gpar(col = "#008a0e", fill = "#008a0e")
)
connect9 <- connectGrob(
  start = pred2,
  end = concl,
  "vertical",
  arrow_obj = custom_arrow
)

# Printing connections
## Note that end-points printed in the console of an IDE does not necessarily 
## reflect their positions in the final printed pdf file
connect1; connect2; connect3; grid.draw(new_connection); 
connect5; connect6; connect7; connect8; connect9
no; yes

flow_chart <- grid.grab()

# Prediction results plot ------------------------------------------------------
# Setting factor levels
res[, ":="(
  Host = factor(Host, levels = c("Lotus", "Hordeum")),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Obs = factor(Obs, levels = c("NPK", "PK", "UF")),
  Pred = factor(Pred, levels = c("UF", "PK", "NPK")),
  Genotype = factor(
    Genotype,
    levels = c("WT", "symrk", "ccamk", "nsp1", "nsp2")
  )
)]
res[, Prediction := Obs == Pred]
genotype_labels_legend <- c(
  "WT" = "WT",
  "symrk" = "*symrk*",
  "ccamk" = "*ccamk*",
  "nsp1" = "*nsp1*",
  "nsp2" = "*nsp2*"
)

bubble_plot <- ggplot(data = res) +
  geom_count(aes(x = Obs, y = Pred), color = "lightgrey") +
  scale_size_continuous(range = c(1.5, 15)) +
  geom_jitter(
    data = res[Prediction == TRUE],
    mapping = aes(x = Obs, y = Pred, fill = Genotype),
    position = position_jitter(width = 0.35, height = 0.35, seed = 1),
    shape = 21,
    stroke = 0.25
  ) +
  geom_jitter(
    data = res[Prediction == FALSE],
    mapping = aes(x = Obs, y = Pred, fill = Genotype),
    position = position_jitter(width = 0.15, height = 0.15, seed = 1),
    shape = 21,
    stroke = 0.25
  ) +
  scale_shape_manual(values = c(21, 21)) +
  theme_bw() +
  labs(x = "Observed", y = "Predicted") +
  guides(size = "none", fill = guide_legend(override.aes = list(size = 3))) +
  facet_grid(Host ~ Compartment) +
  scale_fill_manual(values = colors_genotype, labels = genotype_labels_legend) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -8),
    legend.box.margin = margin(0, 0, 0, -15),
    strip.background = element_rect(colour = NA),
    axis.title.y = element_text(size = 6, family = "Helvetica"),
    axis.title.x = element_text(size = 6, family = "Helvetica"),
    axis.text.y = element_text(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    legend.text = element_markdown(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    legend.title = element_text(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    strip.text = element_text(size = 6, family = "Helvetica", face = "bold"),
    legend.key.size = unit(5, "mm")
  ) +
  NULL

# Bar plot showing taxonomic composition ---------------------------------------
# Setting factor levels
tax_summary[, ":="(
  Order = droplevels(factor(Order, levels = colors_order$Order)),
  Plant = factor(Plant, levels = c("Lotus", "Hordeum")),
  Compartment = factor(Compartment, levels = c("Rhizosphere", "Root")),
  Soil = factor(Soil, levels = c("NPK", "PK", "UF"))
)]
tax_summary <- tax_summary[order(Order)]

bar_plot <- ggplot(data = tax_summary, aes(x = Soil, y = RA, fill = Order)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0.1) +
  facet_wrap2(vars(Plant, Compartment), strip = strip_nested(), nrow = 1) +
  scale_fill_manual(
    values = colors_order$Color,
    breaks = colors_order$Order,
    name = "Bacterial orders"
  ) +
  labs(x = NULL, y = "Cumulative Mean Relative Abundance") +
  guides(fill = guide_legend(nrow = 3, title.position = "top")) +
  scale_y_continuous(expand = c(.001, .001), limits = c(0, 0.38)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, -15),
    legend.margin = margin(t = -5),
    strip.background = element_rect(colour = NA),
    legend.key.size = unit(0.25, 'cm'),
    legend.key.spacing.y = unit(0, 'cm'),
    legend.justification = c(0.75, 0),
    axis.title.y = element_text(size = 6, family = "Helvetica"),
    axis.text.y = element_text(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 6,
      family = "Helvetica",
      colour = "black"
    ),
    legend.text = element_text(size = 6, family = "Helvetica"),
    legend.title = element_text(size = 6, family = "Helvetica"),
    strip.text = element_text(size = 6, family = "Helvetica", face = "bold")
  ) +
  NULL

#  Combining plots -------------------------------------------------------------
accuracy_table <- tableGrob(
  ratio_amount,
  theme = ttheme_default(base_size = 6),
  rows = NULL
)

# First column of whole figure
col1 <- plot_grid(
  flow_chart,
  accuracy_table,
  labels = c("A", "C"),
  rel_heights = c(0.8, 0.2),
  ncol = 1,
  label_size = 12
)

# Second column of whole figure
col2 <- plot_grid(
  bubble_plot,
  bar_plot,
  labels = c("B", "D"),
  rel_heights = c(0.5, 0.5),
  ncol = 1,
  label_size = 12
)

# Combining columns
final_figure <- plot_grid(col1, col2, ncol = 2, label_size = 12)

# Saving plot
ggsave(
  filename = "Figure3_Askov_prediction.pdf",
  plot = final_figure,
  width = 18,
  height = 20,
  units = "cm"
)
