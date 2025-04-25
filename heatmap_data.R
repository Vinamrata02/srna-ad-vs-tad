library(pheatmap)
library(RColorBrewer)

labels <- read_csv("/N/u/vinshar/Quartz/Downloads/ROSMAP_Clinical_Subtypes(in)(1).csv")

# order
labels <- labels[order(labels$Subtype), ]

# get indices
subtypes <- table(labels$Subtype)
indices <- as.vector(cumsum(subtypes))

# for subtype
features <- c("Subtype")
subset <- labels[, features]
subset$Subtype <- as.numeric(as.factor(subset$Subtype))
colnames(subset) <- c("")
heatmap_matrix <- as.matrix(subset)
pheatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = brewer.pal(4, "BrBG"),
         labels_col = "Subtype",
         legend_breaks = c(1, 2, 3, 4),
         legend_labels = c("AsymAD", "LowNFT", "TAD", "CN"),
         gaps_col = indices,
         show_colnames = FALSE,
         angle_col = 0, cellheight = 40, cellwidth = 0.9)


# for braak stage
features <- c("braaksc")
subset <- labels[, features]
subset$braaksc <- as.numeric(as.factor(subset$braaksc))
colnames(subset) <- c("")
heatmap_matrix <- as.matrix(subset)
pheatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = brewer.pal(6, "Dark2"), 
         gaps_col = indices,
         show_colnames = FALSE,
         legend_breaks = c(1, 2, 3, 4, 5, 6), angle_col = 0, cellheight = 40, cellwidth = 0.9)


# for sex
features <- c("sex")
subset <- labels[, features]
subset$sex <- as.numeric(as.factor(subset$sex))
colnames(subset) <- c("")
heatmap_matrix <- as.matrix(subset)
pheatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = c("blue", "green"),
         gaps_col = indices,
         show_colnames = FALSE,
         legend_breaks = c(1, 2),
         legend_labels = c("M", "F"),
         angle_col = 0, cellheight = 40, cellwidth = 0.9)

# for age at death
features <- c("age_death")
subset <- labels[, features]
subset <- subset %>% mutate(age_death = gsub("\\+", "", age_death))
subset$age_death <- as.numeric(subset$age_death)
colnames(subset) <- c("")
heatmap_matrix <- as.matrix(subset)
pheatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("yellow", "purple"))(50),
         gaps_col = indices,
         show_colnames = FALSE,
         angle_col = 0, cellheight = 40, cellwidth = 0.9)

# post_mortem_interval
features <- c("pmi")
subset <- labels[, features]
subset$pmi <- as.numeric(subset$pmi)
colnames(subset) <- c("")
heatmap_matrix <- as.matrix(subset)
pheatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("cyan", "purple"))(10),
         gaps_col = indices,
         show_colnames = FALSE,
         angle_col = 0, cellheight = 40, cellwidth = 0.9)

# congitive differential and neuratic plaque
features <- c("cogdx", "ceradsc")
subset <- labels[, features]
colnames(subset) <- c("", "")
heatmap_matrix <- as.matrix(subset)
pheatmap(t(heatmap_matrix), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("orange", "yellow"))(10),
         gaps_col = indices,
         show_colnames = FALSE,
         angle_col = 0, cellheight = 40, cellwidth = 0.9)

