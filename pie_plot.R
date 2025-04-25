library(Seurat)
library(ggplot2)
library(dplyr)

# extract metadata
#metadata <- object@meta.data

# read labels
labels <- read_csv("/N/u/vinshar/Quartz/Downloads/ROSMAP_Clinical_Subtypes(in)(1).csv")

# join
#metadata <- left_join(metadata, labels)

# remove NA rows
#metadata <- na.omit(metadata)

# pie chart
create_pie_chart <- function(data, category, title) {
  df <- data %>%
    group_by(!!sym(category)) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(100 * count / sum(count), 1))
  
  ggplot(df, aes(x = "", y = count, fill = as.factor(!!sym(category)))) +
    geom_bar(stat = "identity", width = 1) + 
    coord_polar("y", start = 0) + 
    theme_void() +
    theme(legend.position = "right") + 
    labs(fill = category, title = title) + 
    geom_text(aes(label = paste0(count, " (", percentage, "%)")), 
              position = position_stack(vjust = 0.5))
  
}

plot1 <- create_pie_chart(labels, "Study", "Study")
plot2 <- create_pie_chart(labels, "sex", "Sex")
plot3 <- create_pie_chart(labels, "Subtype", "Subtype")
plot4 <- create_pie_chart(labels, "braaksc", "Braak Stage")