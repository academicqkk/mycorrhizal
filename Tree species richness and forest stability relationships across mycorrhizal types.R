# Load required libraries for data manipulation, modeling, and visualization
library(ggplot2)   # For data visualization
library(dplyr)     # For data wrangling
library(forcats)   # For factor level reordering
library(purrr)     # For functional programming tools (e.g., map)
library(broom)     # For converting model outputs to tidy format
library(tidyr)     # For data tidying (e.g., unnesting)

### Load the dataset Mixed_mycorrhizal_associations_and_ectomycorrhizal_dominance ###

# Fit a linear model (Rs_log ~ species_richness_log) for each mycorrhizal type
data_stats_p1 <- data_all %>%
  group_by(Mycorrhizal.type) %>%
  summarise(
    model = list(lm(Rs_log ~ species_richness_log, data = cur_data())),  # Fit model per group
    .groups = 'drop'
  ) %>%
  mutate(
    summary = map(model, broom::glance),        # Extract model summary statistics (R², F-stat, p-value)
    coefficients = map(model, broom::tidy),     # Extract model coefficients (slope, intercept, etc.)
    df_num = map_dbl(model, ~length(coef(.x)) - 1),  # Degrees of freedom numerator (number of predictors)
    df_den = map_dbl(model, ~df.residual(.x))        # Degrees of freedom denominator (residuals)
  ) %>%
  unnest(cols = c(summary, coefficients), names_sep = "_") %>%  # Flatten nested data frames
  filter(coefficients_term == "species_richness_log") %>%        # Keep only the slope term
  select(Mycorrhizal.type, summary_r.squared, summary_statistic, summary_p.value, coefficients_estimate, df_num, df_den) %>%
  rename(  # Rename for clarity
    F_value = summary_statistic, 
    slope = coefficients_estimate, 
    p_value = summary_p.value, 
    r_squared = summary_r.squared
  )


# Create a faceted scatter plot showing the resistance vs. species richness relationship
p1 = ggplot(data_all, aes(
  x = log2(species_richness), 
  y = Rs_log,
  color = Mycorrhizal.type,
  fill = Mycorrhizal.type,
  shape = Mycorrhizal.type
)) +
  # Customize x-axis: set breaks and labels corresponding to richness levels
  scale_x_continuous(
    breaks =  c(0, 1, 2, 3, 4.17), 
    limits = c(0, 4.17), 
    expand = c(0.03, 0.03),
    labels = c("1", "2", "4", "8", "18")
  ) +
  # Customize y-axis
  scale_y_continuous(
    breaks = seq(0, 16, 4), 
    limits = c(0, 16),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  # Facet the plot by mycorrhizal type, in specified order
  facet_wrap(~ fct_relevel(Mycorrhizal.type, "EcM", "AM", "Mixed EcM and AM"), 
             ncol = 1, strip.position = "right") +
  # Customize shape, fill, and color for each group
  scale_shape_manual(name = "Mycorrhizal type", values = c(21, 23, 24)) +
  scale_fill_manual(name = "Mycorrhizal type", 
                    values = c("EcM" = "mediumseagreen",
                               "AM" = "dodgerblue",
                               "Mixed EcM and AM" = "grey60")) + 
  scale_color_manual(name = "Mycorrhizal type", 
                     values = c("EcM" = "mediumseagreen",
                                "AM" = "dodgerblue",
                                "Mixed EcM and AM" = "grey60")) +
  # Add linear regression lines with confidence bands
  geom_smooth(
    aes(group = Mycorrhizal.type), 
    method = 'lm', 
    formula = y ~ x, 
    se = TRUE, 
    size = 0.75, 
    show.legend = FALSE
  ) +
  # Apply a clean theme and customize styling
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 13),
    strip.text.y = element_text(size = 13),
    strip.background.x = element_rect(fill = "#EEEEEF"),
    strip.background.y = element_rect(fill = "#EEEEEF"),
    strip.placement = "outside",
    panel.spacing.x = unit(0.3, "cm"),
    panel.spacing.y = unit(0.3, "cm"),
    plot.margin = unit(c(0.5, 0.2, 0.5, 0.2), "cm"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 13.5, family = "Arial"),
    legend.text = element_text(size = 13.5, family = "Arial"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width  = unit(1.25, 'cm'),
    legend.key = element_rect(fill = "white"),
    legend.background = element_blank()
  ) +
  # Axis labels using mathematical notation
  labs(
    y = expression(Resistance ~ (ln*italic(R)[S])), 
    x = expression(Species ~ richness ~ (italic(S)))
  )
###

### Relationship between species richness and recovery ###

# Fit a linear model (Le_log ~ species_richness_log) for each mycorrhizal type
data_stats_p3 <- data_all %>%
  group_by(Mycorrhizal.type) %>%
  summarise(
    model = list(lm(Le_log ~ species_richness_log, data = cur_data())),
    .groups = 'drop'
  ) %>%
  mutate(
    summary = map(model, broom::glance),
    coefficients = map(model, broom::tidy),
    df_num = map_dbl(model, ~length(coef(.x)) - 1),
    df_den = map_dbl(model, ~df.residual(.x))
  ) %>%
  unnest(cols = c(summary, coefficients), names_sep = "_") %>%
  filter(coefficients_term == "species_richness_log") %>%
  select(Mycorrhizal.type, summary_r.squared, summary_statistic, summary_p.value, coefficients_estimate, df_num, df_den) %>%
  rename(
    F_value = summary_statistic,
    slope = coefficients_estimate,
    p_value = summary_p.value,
    r_squared = summary_r.squared
  )

# Create a faceted scatter plot for recovery vs. species richness
p3 = ggplot(data_all, aes(
  x = log2(species_richness), 
  y = Le_log,
  color = Mycorrhizal.type,
  fill = Mycorrhizal.type,
  shape = Mycorrhizal.type
)) +
  scale_x_continuous(
    breaks =  c(0, 1, 2, 3, 4.17), 
    limits = c(0, 4.17), 
    expand = c(0.03, 0.03),
    labels = c("1", "2", "4", "8", "18")
  ) +
  scale_y_continuous(
    breaks = seq(0, 16, 4), 
    limits = c(0, 16),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  facet_wrap(~ fct_relevel(Mycorrhizal.type, "EcM", "AM", "Mixed EcM and AM"), 
             ncol = 1, strip.position = "right") +
  # Add jitter to reduce point overlap
  geom_jitter(width = 0.10, size = 1, alpha = 0.05) +
  scale_shape_manual(name = "Mycorrhizal type", values = c(21, 23, 24)) +
  scale_fill_manual(name = "Mycorrhizal type", 
                    values = c("EcM" = "mediumseagreen",
                               "AM" = "dodgerblue",
                               "Mixed EcM and AM" = "grey60")) + 
  scale_color_manual(name = "Mycorrhizal type", 
                     values = c("EcM" = "mediumseagreen",
                                "AM" = "dodgerblue",
                                "Mixed EcM and AM" = "grey60")) +
  geom_smooth(
    aes(group = Mycorrhizal.type), 
    method = 'lm', 
    formula = y ~ x, 
    se = TRUE, 
    size = 0.75, 
    show.legend = FALSE
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 13),
    strip.text.y = element_text(size = 13),
    strip.background.x = element_rect(fill = "#EEEEEF"),
    strip.background.y = element_rect(fill = "#EEEEEF"),
    strip.placement = "outside",
    panel.spacing.x = unit(0.3, "cm"),
    panel.spacing.y = unit(0.3, "cm"),
    plot.margin = unit(c(0.5, 0.2, 0.5, 0.2), "cm"),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14, angle = 90, hjust = 0.5),
    plot.title = element_text(size = 18, hjust = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 13.5, family = "Arial"),
    legend.text = element_text(size = 13.5, family = "Arial"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width  = unit(1.25, 'cm'),
    legend.key = element_rect(fill = "white"),
    legend.background = element_blank()
  ) +
  labs(
    y = expression(Recovery ~ (ln*italic(R)[c])), 
    x = expression(Species ~ richness ~ (italic(S)))
  )
###



# Read the dataset speciesrichnessgroups #
### ANALYSIS: Species Richness Groups - RESISTANCE ###
# Calculate summary statistics (mean, lower bound, upper bound) of resistance (Rs_log)
# grouped by species richness and mycorrhizal type
data_summary_p2 <- data_b %>%
  group_by(species_richness, Mycorrhizal.type) %>%
  summarise(mean_Rs_log = mean(Rs_log, na.rm = TRUE),
            lower = mean(Rs_log, na.rm = TRUE) - sd(Rs_log, na.rm = TRUE),
            upper = mean(Rs_log, na.rm = TRUE) + sd(Rs_log, na.rm = TRUE),
            .groups = 'drop')

# Define a function to conduct one-way ANOVA and Tukey HSD for each species richness level
# and return the significance letters based on group comparisons

get_anova_letters_ordered_Rs <- function(data, richness_level) {
  sub_data <- data %>% filter(species_richness == richness_level)
  
  if (length(unique(sub_data$Mycorrhizal.type)) < 2) {
    return(data.frame(
      species_richness = richness_level,
      Mycorrhizal.type = unique(sub_data$Mycorrhizal.type),
      mean_Rs_log = mean(sub_data$Rs_log),
      significance = NA
    ))
  }
  
  # Fit ANOVA
  aov_model <- aov(Rs_log ~ Mycorrhizal.type, data = sub_data)
  tukey <- TukeyHSD(aov_model)
  tukey_df <- as.data.frame(tukey$Mycorrhizal.type)
  tukey_df$Comparison <- rownames(tukey_df)
  pvals <- tukey_df$`p adj`
  names(pvals) <- tukey_df$Comparison
  
  # Get letter groups
  raw_letters <- multcompLetters(pvals)$Letters
  
  # Compute group means
  means <- sub_data %>%
    group_by(Mycorrhizal.type) %>%
    summarise(mean_Rs_log = mean(Rs_log), .groups = "drop")
  
  # Merge letters with means
  result <- means %>%
    mutate(letter_raw = raw_letters[Mycorrhizal.type]) %>%
    arrange(mean_Rs_log)  # Sort by mean ascending
  
  # Re-assign letters: smallest mean gets "a", next gets "b", etc.
  # Note: keep unique raw letter patterns (e.g., "ab") unchanged, but re-map to follow order
  unique_patterns <- unique(result$letter_raw)
  ordered_letters <- letters[1:length(unique_patterns)]
  names(ordered_letters) <- unique_patterns
  
  result <- result %>%
    mutate(
      significance = ordered_letters[letter_raw],
      species_richness = richness_level
    ) %>%
    select(species_richness, Mycorrhizal.type, significance)
  
  return(result)
}

# Apply the ANOVA + Tukey test function across richness levels 1 through 8
results_ordered_Rs <- map_dfr(1:8, ~get_anova_letters_ordered_Rs(data_b, .x))

# Display the resulting significance groupings
print(results_ordered_Rs)

# Set factor levels to ensure consistent ordering in plots
data_b$Mycorrhizal.type <- factor(data_b$Mycorrhizal.type, levels = c("EcM", "AM", "Mixed EcM and AM"))
data_summary_p2$Mycorrhizal.type <- factor(data_summary_p2$Mycorrhizal.type, levels = c("EcM", "AM", "Mixed EcM and AM"))

# Plot violin + error bars + mean points for resistance
p2 = ggplot(data_b, aes(x = factor(species_richness), y = Rs_log)) +
  geom_violin(aes(fill = Mycorrhizal.type), alpha = 0.4, color = NA,
              position = position_dodge(width = 0.8), show.legend = FALSE) +
  
  geom_errorbar(data = data_summary_p2,
                aes(x = factor(species_richness), ymin = lower, ymax = upper,
                    group = Mycorrhizal.type, color = Mycorrhizal.type),
                width = 0, size = 1, position = position_dodge(width = 0.8),
                inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_point(data = data_summary_p2,
             aes(x = factor(species_richness), y = mean_Rs_log,
                 shape = Mycorrhizal.type, color = Mycorrhizal.type, fill = Mycorrhizal.type),
             size = 4, position = position_dodge(width = 0.8), inherit.aes = FALSE) +
  
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, 4),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
  # Manually define colors and shapes for consistency
  scale_fill_manual(name = "Mycorrhizal type", values = c("EcM" = "mediumseagreen",
                                                          "AM" = "dodgerblue",
                                                          "Mixed EcM and AM" = "grey60")) +
  scale_color_manual(name = "Mycorrhizal type", values = c("EcM" = "mediumseagreen",
                                                           "AM" = "dodgerblue",
                                                           "Mixed EcM and AM" = "grey60")) +
  scale_shape_manual(name = "Mycorrhizal type", values = c("EcM" = 21,
                                                           "AM" = 23,
                                                           "Mixed EcM and AM" = 24)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  
  theme_minimal() + theme_bw() +   # Combine clean themes
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(hjust=0.5, size = 14, angle=90),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "bottom",
    legend.title = element_text(size = 13.5, family = "Arial"),
    legend.text = element_text(size = 13.5, family = "Arial"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(1.25, 'cm'),
    legend.background = element_blank()
  ) +
  labs(y = expression(Resistance ~ (ln*italic(R)[S])), 
       x = expression(Species ~ richness ~ (italic(S))))


### ANALYSIS: Species Richness Groups - RECOVERY ###

# Summarize recovery (Le_log) data by richness and mycorrhizal type
data_summary_p4 <- data_b %>%
  group_by(species_richness, Mycorrhizal.type) %>%
  summarise(mean_Le_log = mean(Le_log, na.rm = TRUE),
            lower = mean(Le_log, na.rm = TRUE) - sd(Le_log, na.rm = TRUE),
            upper = mean(Le_log, na.rm = TRUE) + sd(Le_log, na.rm = TRUE),
            .groups = 'drop')

# Define similar ANOVA + Tukey function for Le_log (recovery)
get_anova_letters_ordered_Le <- function(data, richness_level) {
  sub_data <- data %>% filter(species_richness == richness_level)
  
  if (length(unique(sub_data$Mycorrhizal.type)) < 2) {
    return(data.frame(
      species_richness = richness_level,
      Mycorrhizal.type = unique(sub_data$Mycorrhizal.type),
      mean_Le_log = mean(sub_data$Le_log),
      significance = NA
    ))
  }
  
  # Fit ANOVA
  aov_model <- aov(Le_log ~ Mycorrhizal.type, data = sub_data)
  tukey <- TukeyHSD(aov_model)
  tukey_df <- as.data.frame(tukey$Mycorrhizal.type)
  tukey_df$Comparison <- rownames(tukey_df)
  pvals <- tukey_df$`p adj`
  names(pvals) <- tukey_df$Comparison
  
  # Get letter groups
  raw_letters <- multcompLetters(pvals)$Letters
  
  # Compute group means
  means <- sub_data %>%
    group_by(Mycorrhizal.type) %>%
    summarise(mean_Le_log = mean(Le_log), .groups = "drop")
  
  # Merge letters with means
  result_Le <- means %>%
    mutate(letter_raw = raw_letters[Mycorrhizal.type]) %>%
    arrange(mean_Le_log)  # Sort by mean ascending
  
  # Re-assign letters: smallest mean gets "a", next gets "b", etc.
  # Note: keep unique raw letter patterns (e.g., "ab") unchanged, but re-map to follow order
  unique_patterns <- unique(result_Le$letter_raw)
  ordered_letters <- letters[1:length(unique_patterns)]
  names(ordered_letters) <- unique_patterns
  
  result_Le <- result_Le %>%
    mutate(
      significance = ordered_letters[letter_raw],
      species_richness = richness_level
    ) %>%
    select(species_richness, Mycorrhizal.type, significance)
  
  return(result_Le)
}

# Apply across richness levels 1–8
results_ordered_Le <- map_dfr(1:8, ~get_anova_letters_ordered_Le(data_b, .x))

# View significance results for recovery
print(results_ordered_Le)

# Ensure consistent factor levels for plotting
data_b$Mycorrhizal.type <- factor(data_b$Mycorrhizal.type, levels = c("EcM", "AM", "Mixed EcM and AM"))
data_summary_p4$Mycorrhizal.type <- factor(data_summary_p4$Mycorrhizal.type, levels = c("EcM", "AM", "Mixed EcM and AM"))

# Plot violin + error bar + mean point plot for recovery
p4 = ggplot(data_b, aes(x = factor(species_richness), y = Le_log)) +
  geom_violin(aes(fill = Mycorrhizal.type), alpha = 0.4, color = NA,
              position = position_dodge(width = 0.8), show.legend = FALSE) +
  
  geom_errorbar(data = data_summary_p4,
                aes(x = factor(species_richness), ymin = lower, ymax = upper,
                    group = Mycorrhizal.type, color = Mycorrhizal.type),
                width = 0, size = 1, position = position_dodge(width = 0.8),
                inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_point(data = data_summary_p4,
             aes(x = factor(species_richness), y = mean_Le_log,
                 shape = Mycorrhizal.type, color = Mycorrhizal.type, fill = Mycorrhizal.type),
             size = 4, position = position_dodge(width = 0.8), inherit.aes = FALSE) +
  
  scale_y_continuous(breaks = seq(0, 16, 4), limits = c(0, 16),
                     expand = expansion(mult = c(0.02, 0.02))) +
  
  scale_fill_manual(name = "Mycorrhizal type", values = c("EcM" = "mediumseagreen",
                                                          "AM" = "dodgerblue",
                                                          "Mixed EcM and AM" = "grey60")) +
  scale_color_manual(name = "Mycorrhizal type", values = c("EcM" = "mediumseagreen",
                                                           "AM" = "dodgerblue",
                                                           "Mixed EcM and AM" = "grey60")) +
  scale_shape_manual(name = "Mycorrhizal type", values = c("EcM" = 21,
                                                           "AM" = 23,
                                                           "Mixed EcM and AM" = 24)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  
  theme_minimal() + theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(hjust=0.5, size = 14, angle=90),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.position = "bottom",
    legend.title = element_text(size = 13.5, family = "Arial"),
    legend.text = element_text(size = 13.5, family = "Arial"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(1.25, 'cm'),
    legend.background = element_blank()
  ) +
  labs(y = expression(Recovery ~ (ln*italic(R)[c])), 
       x = expression(Species ~ richiness ~ (italic(S))))
