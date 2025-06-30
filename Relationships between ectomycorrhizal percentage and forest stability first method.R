rm(list=ls())
# Load required libraries
library(nlme)       # For linear mixed-effects models
library(purrr)      # For functional programming tools like map functions
library(dplyr)      # For data manipulation
library(sjPlot)     # For plotting mixed models
library(ggplot2)    # For data visualization
library(ggpubr)     # For enhanced ggplot2 functionalities
library(cowplot)    # For combining multiple plots
library(ggnewscale) # For multiple color scales in ggplot
library(tibble)     # For tidy data frames
library(scales)     # For scaling functions in visualizations

### Load the dataset Mixed_mycorrhizal_associations_and_ectomycorrhizal_dominance ###

################### Fit model for Recovery ###################
data_all <- data %>% 
  mutate(LON = ifelse(duplicated(LON) & duplicated(LAT), LON + runif(n(), -0.0001, 0.0001), LON),
         LAT = ifelse(duplicated(LON) & duplicated(LAT), LAT + runif(n(), -0.0001, 0.0001), LAT))

# Fit a linear mixed-effects model to predict log-transformed recovery (Le_log)
# using standardized predictors, including species richness and soil/environmental variables
# Random effects: nested structure by ECOSUBCD and plot_ID
# Spatial autocorrelation structure included via longitude and latitude
model_data_Le <- lme(
  Le_log ~ scale(species_richness_log) + scale(EcM_percentage_square) + scale(EcM_percentage) +
    scale(STDAGE) + scale(SLOPE) + scale(ELEV) + 
    scale(ph) + scale(totaln) + scale(bulk) + scale(clay) +
    scale(pre) + scale(tem) + scale(Rs_log),
  random = ~ 1 | ECOSUBCD/plot_ID,
  correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
  data = data_all
)

# Display model summary and extract fixed effects
summary(model_data_Le)
fixed_effects_all_Le <- fixef(model_data_Le)

# Extract coefficients of the quadratic and linear terms of EcM percentage
a_total <- as.numeric(fixed_effects_all_Le["scale(EcM_percentage_square)"])
b_total <- as.numeric(fixed_effects_all_Le["scale(EcM_percentage)"])

# Calculate the critical point (maximum) of the quadratic curve (in standardized space)
x_critical_centered_total <- -b_total / (2 * a_total)
x_critical_centered_total

# Retrieve mean and standard deviation of original (unscaled) EcM percentage
mean_EcM_total <- mean(data_all$EcM_percentage, na.rm = TRUE)
sd_EcM_total <- sd(data_all$EcM_percentage, na.rm = TRUE)

# Convert the centered critical point back to the original scale
EcM_critical_original_total <- x_critical_centered_total * sd_EcM_total + mean_EcM_total
EcM_critical_original_total

############################################

# Define a function to compute the EcM critical value for a given species richness level
calculate_EcM_critical <- function(species_richness_value) {
  # Subset data for the specified species richness
  data_filtered <- data_all %>% filter(species_richness == species_richness_value)
  
  # Fit a linear mixed-effects model predicting Rl_log
  model_data <- lme(
    Rl_log ~ scale(EcM_percentage_square) + scale(EcM_percentage) +
      scale(STDAGE) + scale(SLOPE) + scale(ELEV) + 
      scale(ph) + scale(totaln) + scale(bulk) + scale(clay) +
      scale(pre) + scale(tem) + scale(Rs_log),
    random = ~ 1 | ECOSUBCD/plot_ID,
    correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
    data = data_filtered
  )
  
  # Extract fixed effect coefficients
  fixed_effects <- fixef(model_data)
  a <- as.numeric(fixed_effects["scale(EcM_percentage_square)"])
  b <- as.numeric(fixed_effects["scale(EcM_percentage)"])
  
  # Compute critical point in scaled space and convert to original scale
  x_critical_centered <- -b / (2 * a)
  mean_EcM <- mean(data_filtered$EcM_percentage, na.rm = TRUE)
  sd_EcM <- sd(data_filtered$EcM_percentage, na.rm = TRUE)
  EcM_critical_original <- x_critical_centered * sd_EcM + mean_EcM
  
  # Return a data frame with species richness and critical EcM value
  return(data.frame(species_richness = species_richness_value, EcM_critical_original = EcM_critical_original))
}

# Compute critical EcM percentage for species richness from 2 to 18
results <- map_dfr(2:18, calculate_EcM_critical)
print(results)

# Compute average critical EcM value across all richness levels
mean_2 = mean(results$EcM_critical_original)

# Create a lollipop plot to visualize the relationship
p_2 = ggplot(results, aes(x = species_richness, y = EcM_critical_original)) +
  geom_hline(yintercept = mean_2, color = "#ad39ad", size = 0.5, linetype = "longdash") +
  geom_segment(aes(xend = species_richness, y = 0, yend = EcM_critical_original),
               color = "grey50", alpha = 0.6, linetype = "dotted", size = 0.8) +
  geom_point(aes(fill = EcM_critical_original), shape = 24, size = 4) +
  scale_fill_gradient2(
    name = "Ectomycorrhizal percentage",
    low = "dodgerblue", mid = "#EBF4F8", high = "mediumseagreen", midpoint = 49.51148
  ) +
  labs(x = expression(Species~richness~(italic(S))),
       y = "The highest recovery \ncorresponds to EcM percentage (%)") +
  scale_x_continuous(breaks = seq(2, 18, 4), limits = c(2, 18),
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0.5, 0.55, 0.5, 0.55), "cm"),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 16, family = "Arial", color = "black", angle = 90),
    axis.title.x = element_text(size = 16, family = "Arial"),
    axis.title.y = element_text(size = 16, family = "Arial"),
    axis.ticks.length = unit(1, "mm"),
    panel.grid = element_blank()
  )

###Fit model for Resistance###

# Fit a linear mixed-effects model predicting Rs_log
model_data_Rs <- lme(
  Rs_log ~ scale(species_richness_log) + scale(EcM_percentage_square) + scale(EcM_percentage) +
    scale(STDAGE) + scale(SLOPE) + scale(ELEV) +
    scale(ph) + scale(totaln) + scale(bulk) + scale(clay) +
    scale(pre) + scale(tem),
  random = ~ 1 | ECOSUBCD/plot_ID,
  correlation = corExp(form = ~ LON + LAT, nugget = TRUE),
  data = data_all
)
summary(model_data_Rs)
fixed_effects_all_Rs <- fixef(model_data_Rs)

# Extract coefficients and compute critical EcM percentage
a_total <- as.numeric(fixed_effects_all_Rs["scale(EcM_percentage_square)"])
b_total <- as.numeric(fixed_effects_all_Rs["scale(EcM_percentage)"])
x_critical_centered_total <- -b_total / (2 * a_total)
mean_EcM_total <- mean(data_all$EcM_percentage, na.rm = TRUE)
sd_EcM_total <- sd(data_all$EcM_percentage, na.rm = TRUE)
EcM_critical_original_total <- x_critical_centered_total * sd_EcM_total + mean_EcM_total
EcM_critical_original_total

# Repeat critical EcM calculation for Rs model by species richness
calculate_EcM_critical <- function(species_richness_value) {
  data_filtered <- data_all %>% filter(species_richness == species_richness_value)
  
  model_data <- lme(
    Rs_log ~ scale(EcM_percentage_square) + scale(EcM_percentage) +
      scale(STDAGE) + scale(SLOPE) + scale(ELEV) + 
      scale(ph) + scale(totaln) + scale(bulk) + scale(clay) +
      scale(pre) + scale(tem),
    random = ~ 1 | ECOSUBCD/plot_ID,
    correlation = corExp(form = ~ LON + LAT),
    data = data_filtered
  )
  
  fixed_effects <- fixef(model_data)
  a <- as.numeric(fixed_effects["scale(EcM_percentage_square)"])
  b <- as.numeric(fixed_effects["scale(EcM_percentage)"])
  x_critical_centered <- -b / (2 * a)
  mean_EcM <- mean(data_filtered$EcM_percentage, na.rm = TRUE)
  sd_EcM <- sd(data_filtered$EcM_percentage, na.rm = TRUE)
  EcM_critical_original <- x_critical_centered * sd_EcM + mean_EcM
  
  return(data.frame(species_richness = species_richness_value, EcM_critical_original = EcM_critical_original))
}

# Compute critical EcM percentages for Rs across species richness values
results_2 <- map_dfr(2:18, calculate_EcM_critical)
print(results_2)

mean_1 = mean(results_2$EcM_critical_original)

# Create lollipop plot for Rs model
p_1 = ggplot(results_2, aes(x = species_richness, y = EcM_critical_original)) +
  geom_hline(yintercept = mean_1, color = "#ad39ad", size = 0.5, linetype = "longdash") +
  geom_segment(aes(xend = species_richness, y = 0, yend = EcM_critical_original),
               color = "grey50", alpha = 0.6, linetype = "dotted", size = 0.8) +
  geom_point(aes(fill = EcM_critical_original), shape = 24, size = 4) +
  scale_fill_gradient2(
    name = "Ectomycorrhizal percentage",
    low = "dodgerblue", mid = "#EBF4F8", high = "mediumseagreen", midpoint = 49.51148
  ) +
  labs(
    x = expression(Species~richness~(italic(S))),
    y = "The highest resistance\ncorresponds to EcM percentage (%)"
  ) +
  scale_x_continuous(breaks = seq(2, 18, 4), limits = c(2, 18),
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0.5, 0.55, 0.5, 0.55), "cm"),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 16, family = "Arial", color = "black", angle = 90),
    axis.title.x = element_text(size = 16, family = "Arial"),
    axis.title.y = element_text(size = 16, family = "Arial"),
    axis.ticks.length = unit(1, "mm"),
    panel.grid = element_blank()
  )
###

### Load the dataset EcMpercentagemeansgrouped_resistance ###
###EcM percentage groups-resistance analysis###
# Group data by EcM percentage intervals (each 1%)
# Then calculate the group-wise mean and standard error for Rs_log,
# along with mean EcM percentage and species richness
# Ensure numeric type for plotting
grouped_means_Rs$group = as.numeric(grouped_means_Rs$group)
grouped_means_Rs$mean_EcM = as.numeric(grouped_means_Rs$mean_EcM)

# Check species richness range for reference
range(grouped_means_Rs$mean_species_richness)

# Plot Rs against EcM percentage groups with error bars and smooth trend
p1 = ggplot(data_all, aes(x = EcM_percentage)) +
  
  # Add vertical reference line at 80% EcM (optional threshold)
  geom_segment(aes(x = 80, xend = 80, y = 1.5, yend = 5.3), 
               color = "grey60", linetype = "dotted", size = 1) +
  
  # Add error bars (± standard error of Rs_log)
  geom_errorbar(data = grouped_means_Rs, 
                aes(x = group, ymin = mean_Rs_log - se_Rs_log, ymax = mean_Rs_log + se_Rs_log),
                width = 0, color = "lightgrey") +
  
  # Add points showing mean Rs_log per EcM group, colored and sized by group and richness
  geom_point(data = grouped_means_Rs, 
             aes(x = group, y = mean_Rs_log, 
                 fill = group, size = mean_species_richness), 
             alpha = 0.95, shape = 21) +  
  
  # Gradient fill and size scale
  scale_fill_gradient2(name = "EcM percentage (%)",
                       mid = "#EBF4F8", low = "dodgerblue", high = "mediumseagreen",
                       midpoint = 49.51148) +
  scale_size("Species richness", limits = c(1.4, 6.2), range = c(2, 6)) +
  
  # Add a quadratic GAM smoothing curve
  geom_smooth(data = grouped_means_Rs, 
              aes(x = group, y = mean_Rs_log), 
              method = "gam", formula = y ~ poly(x, 2), 
              color = "black", fill = "lightgrey", alpha = 0.6, size = 0.6, se = TRUE) + 
  
  # X and Y axis formatting
  scale_x_continuous(breaks = seq(0, 100, 10), expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(breaks = c(2, 3, 4, 5, 6, 7), limits = c(1.5, 7.5), expand = expansion(mult = c(0, 0))) +
  
  # Axis and label styling
  labs(x = "EcM percentage (%)", y = expression(Resistance~(italic(R)[s]))) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0.5, 0.55, 0.5, 0.55), "cm"),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 16, family = "Arial", color = "black", angle = 90),
    axis.title = element_text(size = 16, family = "Arial", color = "black"),
    axis.ticks.length = unit(1, "mm"),
    strip.text = element_text(size = 16, family = "Arial", color = "black"),
    legend.spacing.x = unit(2, "cm"),
    legend.spacing.y = unit(0.5, "cm")
  ) +
  # Customize legends
  guides(
    size = guide_legend(title.position = "top", title.hjust = 0.5, label.position = "bottom", keywidth = unit(0.8, "cm")),
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5, label.position = "bottom", barwidth = unit(6.5, "cm"))
  )
###

### Load the dataset EcMpercentagemeansgrouped_recovery ###
###EcM percentage groups-recovery analysis###
# Group data by EcM percentage bins (same as above)
# Compute group-wise mean and standard error for Le_log and associated metadata
# Convert group and mean EcM percentage to numeric for plotting
grouped_means_Le$group = as.numeric(grouped_means_Le$group)
grouped_means_Le$mean_EcM = as.numeric(grouped_means_Le$mean_EcM)
# Plot Le against EcM percentage groups with smoothing and styling
p2 = ggplot(data_all, aes(x = EcM_percentage)) +
  
  # Reference line at 81% EcM (could indicate a threshold)
  geom_segment(aes(x = 81, xend = 81, y = 1.5, yend = 5.7), 
               color = "grey60", linetype = "dotted", size = 1) +
  
  # Error bars for Le_log ± SE
  geom_errorbar(data = grouped_means_Le, 
                aes(x = group, ymin = mean_Le_log - se_Le_log, ymax = mean_Le_log + se_Le_log),
                width = 0, color = "lightgrey") +
  
  # Plot group means with fill and size encoding
  geom_point(data = grouped_means_Le, 
             aes(x = group, y = mean_Le_log, fill = group, size = mean_species_richness),
             alpha = 0.95, shape = 21) +
  
  # Color and size scales
  scale_fill_gradient2(name = "Ectomycorrhizal percentage",
                       mid = "#EBF4F8", low = "dodgerblue", high = "mediumseagreen",
                       midpoint = 49.51148) +
  scale_size("Species richness", limits = c(1.4, 6.2), range = c(2, 6)) +
  
  # Quadratic smoothing of mean Le_log
  geom_smooth(data = grouped_means_Le, 
              aes(x = group, y = mean_Le_log),
              method = "gam", formula = y ~ poly(x, 2),
              color = "black", fill = "lightgrey", alpha = 0.6, size = 0.6, se = TRUE) +
  
  # Axes
  scale_x_continuous(breaks = seq(0, 100, 10), expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(breaks = seq(2, 8, 2), limits = c(1.5, 8.5), expand = expansion(mult = c(0, 0))) +
  
  # Labels and theme
  labs(x = "EcM percentage (%)", y = expression(Recovery~(italic(R)[c]))) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0.5, 0.55, 0.5, 0.55), "cm"),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
    axis.text.y = element_text(size = 16, family = "Arial", color = "black", angle = 90),
    axis.title = element_text(size = 16, family = "Arial", color = "black"),
    axis.ticks.length = unit(1, "mm"),
    strip.text = element_text(size = 16, family = "Arial", color = "black"),
    legend.spacing.x = unit(2, "cm"),
    legend.spacing.y = unit(0.5, "cm")
  ) +
  guides(
    size = guide_legend(title.position = "top", title.hjust = 0.5, label.position = "bottom", keywidth = unit(0.8, "cm")),
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5, label.position = "bottom", barwidth = unit(6.5, "cm"))
  )
###