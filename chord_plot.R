library(circlize)

# assumes deg and deg_count as in memory.. clean_data.R

# cluster names
cluster_names <- names(deg_count)

# create nxn grid
chord_data <- expand.grid(from = cluster_names, to = cluster_names, stringsAsFactors = FALSE)

# assign values based on the originating cluster (from)
chord_data <- chord_data %>% mutate(value = deg_count[from])

# define colors for each cluster
col_fun <- colorRamp2(c(min(deg_count), max(deg_count)), c("blue", "red"))
cluster_colors <- setNames(rainbow(length(cluster_names)), cluster_names)

# Draw chord diagram
chordDiagram(chord_data, grid.col = cluster_colors)

# Add legend
legend("bottomright", legend = cluster_names, col = cluster_colors, pch = 15, cex = 0.7)

# Print the final data frame
print(chord_data)