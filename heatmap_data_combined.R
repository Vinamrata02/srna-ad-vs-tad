library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(tidyverse)
library(grid)

# Load the data
labels <- read_csv("/N/u/vinshar/Quartz/Downloads/ROSMAP_Clinical_Subtypes(in)(1).csv")

# Define a function to create heatmaps with consistent column order
create_heatmap <- function(data, feature, colors, legend_breaks, legend_labels = NULL, title) {
  subset <- labels[, feature, drop = FALSE]
  
  # Convert categorical variables to numeric factors
  if (feature %in% c("Subtype", "sex", "braaksc")) {
    subset[[feature]] <- as.numeric(as.factor(subset[[feature]]))
  } else if (feature == "age_death") {
    subset <- subset %>% mutate(age_death = gsub("\\+", "", age_death))
    subset$age_death <- as.numeric(subset$age_death)
  } else {
    subset[[feature]] <- as.numeric(subset[[feature]])
  }
  
  # Ensure all heatmaps have the same column order
  subset <- subset[order(row.names(subset)), , drop = FALSE]
  
  colnames(subset) <- title
  heatmap_matrix <- as.matrix(subset)
  
  # Return pheatmap object with clustering disabled
  pheatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
           color = colors,
           legend_breaks = legend_breaks,
           legend_labels = legend_labels,
           angle_col = 0, cellheight = 15, silent = TRUE, border_color = NA)
}

# Generate heatmaps ensuring consistent column order
p1 <- create_heatmap(labels, "Subtype", brewer.pal(4, "BrBG"), c(1, 2, 3, 4), c("AsymAD", "LowNFT", "TAD", "CN"), "Subtype")
p2 <- create_heatmap(labels, "braaksc", brewer.pal(6, "Dark2"), c(1, 2, 3, 4, 5, 6), NULL, "Braak Stage")
p3 <- create_heatmap(labels, "sex", c("blue", "green"), c(1, 2), c("male", "female"), "Sex")
p4 <- create_heatmap(labels, "age_death", colorRampPalette(c("yellow", "purple"))(50), NULL, NULL, "Age at Death")
p5 <- create_heatmap(labels, "pmi", colorRampPalette(c("cyan", "purple"))(10), NULL, NULL, "Post Mortem Interval")
p6 <- create_heatmap(labels, "cogdx", colorRampPalette(c("orange", "yellow"))(10), NULL, NULL, "Change in Cognitive Function")
p7 <- create_heatmap(labels, "ceradsc", colorRampPalette(c("orange", "yellow"))(10), NULL, NULL, "Neuritic Plaque Burden")

# Align all heatmaps precisely in a grid
grid.arrange(grobs = list(p1$gtable, p2$gtable, p3$gtable, p4$gtable, p5$gtable, p6$gtable, p7$gtable), ncol = 1, respect = TRUE)