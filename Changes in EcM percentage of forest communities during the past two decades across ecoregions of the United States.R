###
library(ggplot2)   # For data visualization
library(dplyr)     # For data wrangling
library(forcats)   # For factor level reordering
library(purrr)     # For functional programming tools (e.g., map)
library(broom)     # For converting model outputs to tidy format
library(tidyr)     # For data tidying (e.g., unnesting)

# Load INVYR_percentage_change #

# Calculate the mean and standard deviation of percentage change for each ecoregion
eco_change_eco_men <- INVYR_plot_Mycorrhizal %>%
  group_by(ECOSUBCD) %>%
  summarise(
    mean = mean(chage_percentage, na.rm = TRUE),  # Average percentage change
    sd = sd(chage_percentage, na.rm = TRUE)       # Standard deviation of percentage change
  )
print(eco_change_eco_men)

# Calculate the number of plots with increased or decreased EcM (ectomycorrhizal) proportions by ecological subsection
proportions_EcM <- INVYR_plot_Mycorrhizal %>%
  group_by(ECOSUBCD) %>%
  summarise(
    increase = sum(chage_percentage > 0),
    decrease = sum(chage_percentage < 0)
  ) %>%
  
  # Calculate the percentage of plots that show a decrease in EcM
  mutate(decrease_percentage = decrease / (increase + decrease) * 100) %>%
  
  # Classify each ecological subsection into a change group based on the percentage of decreasing plots
  mutate(change_group = case_when(
    decrease_percentage < 40 ~ "<40",
    decrease_percentage >= 40 & decrease_percentage < 50 ~ "40-50",
    decrease_percentage >= 50 & decrease_percentage < 60 ~ "50-60",
    decrease_percentage >= 60 & decrease_percentage < 70 ~ "60-70",
    decrease_percentage >= 70 ~ ">70",
    TRUE ~ NA_character_
  )) %>%
  
  ungroup()

# Perform paired Wilcoxon signed-rank tests and determine significance of change
trend_results <- INVYR_plot_Mycorrhizal %>%
  group_by(ECOSUBCD) %>%
  summarise(
    # Calculate p-value using the Wilcoxon test for paired samples
    p_value = wilcox.test(max_INVYR_percentage, min_INVYR_percentage, paired = TRUE)$p.value,
    
    # Compute the average difference between max and min INVYR percentages
    mean_diff = mean(max_INVYR_percentage - min_INVYR_percentage, na.rm = TRUE),
    
    # Categorize the trend based on statistical significance and direction of change
    trend = case_when(
      p_value < 0.05 & mean_diff > 0 ~ "Significant Increase",
      p_value < 0.05 & mean_diff < 0 ~ "Significant Decrease",
      TRUE ~ "No Significant Change"
    ),
    
    .groups = "drop"
  )

# Display the summarized trend results
print(trend_results)
# Define a custom color mapping for each ECOSUBCD code
color_mapping <- c(
  "242" = "#F1CF70", "333" = "#CA6C21", "332" = "#F1CF70", "251" = "#B4BA7C",
  "212" = "#F1CF70", "211" = "#B4BA7C", "334" = "#F1CF70", "342" = "#F1CF70",
  "331" = "#B4BA7C", "221" = "#F1CF70", "263" = "#963F0C", "261" = "#F1CF70",
  "222" = "#F1CF70", "341" = "#728E53", "262" = "#F1CF70", "322" = "#728E53",
  "223" = "#CA6C21", "234" = "#B4BA7C", "313" = "#F1CF70", "231" = "#F1CF70",
  "255" = "#CA6C21", "321" = "#963F0C", "232" = "#B4BA7C", "315" = "#963F0C",
  "411" = "#B4BA7C"
)

# Assign a fill color to each ecoregion based on the color mapping
trend_results$fill <- color_mapping[as.character(trend_results$ECOSUBCD)]

# Define the desired order for color appearance in the plot
sorted_colors <- c("#728E53", "#B4BA7C", "#F1CF70", "#CA6C21", "#963F0C")

# Convert the fill variable to a factor with the specified order
trend_results$fill <- factor(trend_results$fill, levels = sorted_colors)

# Count the number of ecoregions in each trend category
trend_counts <- trend_results %>%
  group_by(trend) %>%
  summarise(count = n())

# Merge the count data with the original trend results for plotting
data_plot <- trend_results %>%
  left_join(trend_counts, by = "trend")

# Reorder the trend categories for the x-axis
data_plot$trend <- factor(data_plot$trend, 
                          levels = c("Significant Decrease",
                                     "Significant Increase",
                                     "No Significant Change"))

# Create a bar plot showing the number of ecoregions by trend and color
p_bar <- ggplot(data_plot, aes(x = trend, fill = fill, color = fill)) +
  geom_bar(alpha = 0.85, width = 0.5) +
  scale_fill_manual(values = sorted_colors) +
  scale_color_manual(values = sorted_colors) +
  labs(x = NULL, y = "Number of Ecoregions", fill = "Color Order") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  scale_y_continuous(
    breaks = seq(0, 10, 1), limits = c(0, 10),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    labels = c("Significant\ndecrease",
               "Significant\nincrease",
               "No significant\nchange")
  ) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.4),
    legend.position = "none",
    axis.text.x = element_text(hjust = 0.5, vjust = 0.5, angle = 90, size = 16, color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(hjust = 0.5, size = 16, angle = 90, margin = margin(r = 5)),
    strip.text = element_text(size = 16, family = "Arial", color = "black")
  )

# Display the plot
p_bar
