# =========================================================
#
# Dependencies:
#   dplyr, ggplot2, ggpubr (or rstatix for stat_pvalue_manual), grid (unit).
# =========================================================

### Load the dataset Arbuscularmycorrhizalforests, Ectomycorrhizalforests, Mixedmycorrhizalforests ###

data_EcM_AM = filtered_data_EcM_AM %>% filter(species_richness == 1)
data_AM_1 = filtered_data_AM %>% filter(species_richness == 1)
data_EcM_1 = filtered_data_EcM %>% filter(species_richness == 1)

# ---------------------------------------------------------
# Sample size counts per group (for optional annotation)
# ---------------------------------------------------------
plot_num_total <- data.frame(
  type = c("AM","EcM","EcM_AM"),
  lbl  = c(dim(data_AM_1)[1], dim(data_EcM_1)[1], dim(data_EcM_AM_1)[1])
)

# ---------------------------------------------------------
# Significance annotations for Rs
# ---------------------------------------------------------
my_comparisons_Rs <- list(c("AM","EcM"), c("EcM","EcM_AM"), c("AM","EcM_AM"))

my_comparisons_labels <- data.frame(
  group1     = c("AM", "AM", "EcM"),
  group2     = c("EcM", "EcM_AM", "EcM_AM"),
  y.position = c(18, 20, 19),   # vertical position of significance bars
  label      = c("***", "***", "***")  # custom significance stars
)

# Fix factor levels for consistent order on x-axis
data_AM_EcM_EcM_AM$type <- factor(
  data_AM_EcM_EcM_AM$type,
  levels = c("AM","EcM","EcM_AM")
)

# ---------------------------------------------------------
# Summary statistics for Rs_log (mean ± SD)
# ---------------------------------------------------------
data_summary_Rs <- data_AM_EcM_EcM_AM %>%
  group_by(type) %>%
  summarise(
    mean_Rs_log = mean(Rs_log, na.rm = TRUE),
    lower       = mean(Rs_log, na.rm = TRUE) - sd(Rs_log, na.rm = TRUE),
    upper       = mean(Rs_log, na.rm = TRUE) + sd(Rs_log, na.rm = TRUE),
    .groups     = 'drop'
  )

# ---------------------------------------------------------
# Plot A: Rs (Resistance) violin + group means + error bars + significance
# ---------------------------------------------------------
pa = ggplot(data_AM_EcM_EcM_AM, aes(x = type, y = Rs_log)) +
  geom_violin(aes(fill = type), alpha = 0.4, color = NA,
              position = position_dodge(width = 0.75), show.legend = FALSE) +
  
  geom_errorbar(data = data_summary_Rs,
                aes(x = type, ymin = lower, ymax = upper, group = type, color = type),
                width = 0, size = 1, position = position_dodge(width = 0.75),
                inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_point(data = data_summary_Rs,
             aes(x = type, y = mean_Rs_log, shape = type, color = type, fill = type),
             size = 4, position = position_dodge(width = 0.75),
             inherit.aes = FALSE, show.legend = FALSE) +
  
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5)) +
  
  scale_fill_manual(name = "Extrem\nclimate",
                    values = c("dodgerblue","mediumseagreen","grey60"),
                    labels = c("AM", "EcM", "Mixed EcM and AM")) +
  scale_color_manual(name = "Extrem\nclimate",
                     values = c("dodgerblue","mediumseagreen","grey60"),
                     labels = c("AM", "EcM", "Mixed EcM and AM")) +
  scale_shape_manual(name = "Extrem\nclimate",
                     values = c("AM" = 21, "EcM" = 23, "EcM_AM" = 24),
                     labels = c("AM\nmonospecific forests",
                                "EcM\nmonospecific forests",
                                "Mixed EcM and AM\nmonospecific forests")) +
  
  scale_x_discrete(name = NULL,
                   breaks = c("AM","EcM","EcM_AM"),
                   labels = c("AM\nmonospecific forests",
                              "EcM\nmonospecific forests",
                              "Mixed EcM and AM\nmonospecific forests")) +
  
  stat_pvalue_manual(my_comparisons_labels,
                     label = "label",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.01, size = 6) +
  
  theme_minimal() + theme_bw() +
  theme(
    strip.text.x = element_text(size = 13),
    strip.text.y = element_text(size = 13),
    strip.background.x = element_rect(color="#EEEEEF", fill="#EEEEEF"),
    strip.background.y = element_rect(color="#EEEEEF", fill="#EEEEEF"),
    strip.placement = "outside",
    panel.spacing.x = unit(0.3, "cm"),
    panel.spacing.y = unit(0.3, "cm"),
    plot.margin = unit(c(0.5,0.2,0.5,0.2), "cm"),
    
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(hjust = 0.5, size = 14, angle = 90),
    plot.title   = element_text(hjust = 0.5, size = 18),
    
    legend.position = "bottom",
    legend.title    = element_text(hjust = 0.5, size = 13.5, family = "Arial", color = "black"),
    legend.text     = element_text(size = 13.5, family = "Arial", color = "black"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width  = unit(1.25, 'cm'),
    legend.key = element_rect(fill = "white"),
    legend.background = element_blank()
  ) +
  labs(y = expression(Resistance ~ (ln*italic(R)[s])),
       x = NULL)

# ---------------------------------------------------------
# Significance annotations for Le (Recovery)
# ---------------------------------------------------------
my_comparisons_Le <- list(c("AM","EcM"), c("EcM","EcM_AM"), c("AM","EcM_AM"))

my_comparisons_labels_Le <- data.frame(
  group1     = c("AM", "AM", "EcM"),
  group2     = c("EcM", "EcM_AM", "EcM_AM"),
  y.position = c(14.4, 16, 15.2),
  label      = c("***", "***", "**")
)

data_AM_EcM_EcM_AM$type <- factor(
  data_AM_EcM_EcM_AM$type,
  levels = c("AM","EcM","EcM_AM")
)

# ---------------------------------------------------------
# Summary statistics for Le_log (mean ± SD)
# NOTE: in your code, SD is taken from Rs_log instead of Le_log — likely a copy-paste typo.
# Keeping your original line here, but see "Optional fix" below.
# ---------------------------------------------------------
data_summary_Le <- data_AM_EcM_EcM_AM %>%
  group_by(type) %>%
  summarise(
    mean_Le_log = mean(Le_log, na.rm = TRUE),
    lower       = mean(Le_log, na.rm = TRUE) - sd(Le_log, na.rm = TRUE),  # <- likely typo
    upper       = mean(Le_log, na.rm = TRUE) + sd(Le_log, na.rm = TRUE),  # <- likely typo
    .groups     = 'drop'
  )

# ---------------------------------------------------------
# Plot B: Le (Recovery)
# ---------------------------------------------------------
pb = ggplot(data_AM_EcM_EcM_AM, aes(x = type, y = Le_log)) +
  geom_violin(aes(fill = type), alpha = 0.4, color = NA,
              position = position_dodge(width = 0.75), show.legend = FALSE) +
  
  geom_errorbar(data = data_summary_Le,
                aes(x = type, ymin = lower, ymax = upper, group = type, color = type),
                width = 0, size = 1, position = position_dodge(width = 0.75),
                inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_point(data = data_summary_Le,
             aes(x = type, y = mean_Le_log, shape = type, color = type, fill = type),
             size = 4, position = position_dodge(width = 0.75),
             inherit.aes = FALSE, show.legend = FALSE) +
  
  scale_y_continuous(limits = c(0, 16), breaks = seq(0, 16, by = 4)) +
  
  scale_fill_manual(name = "Extrem\nclimate",
                    values = c("dodgerblue","mediumseagreen","grey60"),
                    labels = c("AM", "EcM", "Mixed EcM and AM")) +
  scale_color_manual(name = "Extrem\nclimate",
                     values = c("dodgerblue","mediumseagreen","grey60"),
                     labels = c("AM", "EcM", "Mixed EcM and AM")) +
  scale_shape_manual(name = "Extrem\nclimate",
                     values = c("AM" = 21, "EcM" = 23, "EcM_AM" = 24),
                     labels = c("AM\nmonospecific forests",
                                "EcM\nmonospecific forests",
                                "Mixed EcM and AM\nmonospecific forests")) +
  
  scale_x_discrete(name = NULL,
                   breaks = c("AM","EcM","EcM_AM"),
                   labels = c("AM\nmonospecific forests",
                              "EcM\nmonospecific forests",
                              "Mixed EcM and AM\nmonospecific forests")) +
  
  stat_pvalue_manual(my_comparisons_labels_Le,
                     label = "label",
                     xmin = "group1", xmax = "group2",
                     y.position = "y.position",
                     tip.length = 0.0025, size = 6) +
  
  theme_minimal() + theme_bw() +
  theme(
    strip.text.x = element_text(size = 13),
    strip.text.y = element_text(size = 13),
    strip.background.x = element_rect(color="#EEEEEF", fill="#EEEEEF"),
    strip.background.y = element_rect(color="#EEEEEF", fill="#EEEEEF"),
    strip.placement = "outside",
    panel.spacing.x = unit(0.3, "cm"),
    panel.spacing.y = unit(0.3, "cm"),
    plot.margin = unit(c(0.5,0.2,0.5,0.2), "cm"),
    
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(hjust = 0.5, size = 14, angle = 90),
    plot.title   = element_text(hjust = 0.5, size = 18),
    
    legend.position = "bottom",
    legend.title    = element_text(hjust = 0.5, size = 13.5, family = "Arial", color = "black"),
    legend.text     = element_text(size = 13.5, family = "Arial", color = "black"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width  = unit(1.25, 'cm'),
    legend.key = element_rect(fill = "white"),
    legend.background = element_blank()
  ) +
  labs(y = expression(Recovery ~ (ln*italic(R)[c])),
       x = NULL)