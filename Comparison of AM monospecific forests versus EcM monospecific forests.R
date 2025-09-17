# Load necessary libraries
library(ggrepel)     # For better label placement in ggplot2
library(ggplot2)     # Core plotting package
library(dplyr)       # Data manipulation
library(lubridate)   # Date/time functions (not used in this script)
library(grid)        # For unit and graphical layout control

### Load the dataset monospecificforests_pairedcomparisons ###

### --- RESISTANCE ANALYSIS --- ###
# Extract resistance-related variables and classify differences
data_Rs = pair_data %>% 
  dplyr::select(c('ECOSUBCD', 'AM_Rs', 'AM_Rs_sd', 'EcM_Rs', 'EcM_Rs_sd')) %>%
  mutate(delt_Rs = case_when(
    EcM_Rs > AM_Rs ~ "large_EcM",   # EcM resistance is greater
    EcM_Rs < AM_Rs ~ "small_EcM",   # AM resistance is greater
    TRUE ~ "equal"                 # No difference
  ))

# Count the number of ecoregions in each group
dim(data_Rs[data_Rs$delt_Rs == "large_EcM",])[1]  # Number of regions where EcM > AM
dim(data_Rs[data_Rs$delt_Rs == "small_EcM",])[1]  # Number of regions where AM > EcM

# Calculate percentages for annotations in the plot
Rs_down = paste0(
  format(round(dim(data_Rs[data_Rs$delt_Rs == "large_EcM",])[1] / 
                 (dim(data_Rs[data_Rs$delt_Rs == "large_EcM",])[1] +
                    dim(data_Rs[data_Rs$delt_Rs == "small_EcM",])[1]) * 100, 0), nsmall = 0), "%")
Rs_up = paste0(
  format(round(dim(data_Rs[data_Rs$delt_Rs == "small_EcM",])[1] / 
                 (dim(data_Rs[data_Rs$delt_Rs == "large_EcM",])[1] +
                    dim(data_Rs[data_Rs$delt_Rs == "small_EcM",])[1]) * 100, 0), nsmall = 0), "%")

# Define triangular areas for highlighting performance zones
trsup_1 <- data.frame(x=c(0.5, 0.5, 11.5), y=c(0.5, 11.5, 11.5))  # EcM better
trinf_1 <- data.frame(x=c(0.5, 11.5, 11.5), y=c(0.5, 0.5, 11.5))  # AM better

# Separate datasets based on performance
data_Rs_sin = data_Rs[data_Rs$delt_Rs == "small_EcM",]
data_Rs_com = data_Rs[data_Rs$delt_Rs == "large_EcM",]

# Define color schemes
dotCOLS = c("dodgerblue", "mediumseagreen")
barCOLS = c("dodgerblue", "mediumseagreen")

# Set factor levels for consistent plotting
data_Rs$delt_Rs <- factor(data_Rs$delt_Rs, levels = c("small_EcM", "large_EcM"))

# Merge with additional plot count data
data_Rs2 = data_Rs %>% left_join(eco_plot_num, by = "ECOSUBCD")

# Create the resistance scatterplot
p1_1 = ggplot(data_Rs2, aes(x = EcM_Rs, y = AM_Rs)) +
  geom_errorbar(aes(ymin = AM_Rs - AM_Rs_sd, ymax = AM_Rs + AM_Rs_sd),
                width = 0.1, color = "grey50", size = 0.4, alpha = 0.5) +
  geom_errorbarh(aes(xmin = EcM_Rs - EcM_Rs_sd, xmax = EcM_Rs + EcM_Rs_sd),
                 height = 0.15, color = "grey50", size = 0.4, alpha = 0.5) +
  geom_point(aes(fill = delt_Rs, size = ECO_count), 
             color = "black", shape = 21, alpha = 0.85) +
  scale_fill_manual(values = barCOLS, guide = "none") +
  scale_size_continuous(range = c(3, 12),
                        name = "Number of plots\nin ecoregion",
                        breaks = c(50, 1000, 2000, 3000),
                        guide = guide_legend(direction = "horizontal")) +
  scale_y_continuous(name = "Resistance of\nAM monospecific forests",
                     limits = c(0.5, 11.5), breaks = seq(1, 11, 2),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(name = "Resistance of\nEcM monospecific forests",
                     limits = c(0.5, 11.5), breaks = seq(1, 11, 2),
                     expand = expansion(mult = c(0, 0))) +
  geom_abline(slope = 1, linetype = "dashed", color = "black", size = 0.6, alpha = 0.8) +
  geom_polygon(aes(x = x, y = y), data = trsup_1, fill = "dodgerblue", alpha = 0.1) +
  geom_polygon(aes(x = x, y = y), data = trinf_1, fill = "mediumseagreen", alpha = 0.1) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.3), "cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 15, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width  = unit(1.5, 'cm'),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.6),
        axis.text = element_text(size = 15),
        axis.text.y = element_text(hjust = 0.5, angle = 90),
        axis.ticks.length = unit(1, "mm"),
        strip.text = element_text(size = 15, family = "Arial", color = "black")) +
  annotate(geom = "text", x = 6, y = 1, label = Rs_down, size = 6, color = "mediumseagreen") +
  annotate(geom = "text", x = 6, y = 11, label = Rs_up, size = 6, color = "dodgerblue")

p1_1  # Print resistance plot


### --- RECOVERY ANALYSIS --- ###

# Extract recovery-related variables and classify differences
data_Le = pair_data %>%
  dplyr::select(c('ECOSUBCD', 'AM_Le', 'EcM_Le', 'AM_Le_sd', 'EcM_Le_sd')) %>%
  mutate(delt_Le = case_when(
    EcM_Le > AM_Le ~ "large_EcM",   # EcM has better recovery
    EcM_Le < AM_Le ~ "small_EcM",   # AM has better recovery
    TRUE ~ "equal"                  # Equal recovery
  ))

# Check variable ranges (optional debugging)
max(data_Le$AM_Le); min(data_Le$AM_Le)
max(data_Le$EcM_Le); min(data_Le$EcM_Le)

# Count number of regions in each group
dim(data_Le[data_Le$delt_Le == "large_EcM",])[1]
dim(data_Le[data_Le$delt_Le == "small_EcM",])[1]

# Compute percentages for annotations
Le_up = paste0(
  format(round(dim(data_Le[data_Le$delt_Le == "large_EcM",])[1] /
                 (dim(data_Le[data_Le$delt_Le == "large_EcM",])[1] +
                    dim(data_Le[data_Le$delt_Le == "small_EcM",])[1]) * 100, 0), nsmall = 0), "%")
Le_down = paste0(
  format(round(dim(data_Le[data_Le$delt_Le == "small_EcM",])[1] /
                 (dim(data_Le[data_Le$delt_Le == "large_EcM",])[1] +
                    dim(data_Le[data_Le$delt_Le == "small_EcM",])[1]) * 100, 0), nsmall = 0), "%")

# Define triangular areas
trsup_2 <- data.frame(x = c(1.5, 1.5, 8.5), y = c(1.5, 8.5, 8.5))
trinf_2 <- data.frame(x = c(1.5, 8.5, 8.5), y = c(1.5, 1.5, 8.5))

# Prepare subset data
data_Le_sin = data_Le[data_Le$delt_Le == "small_EcM",]
data_Le_com = data_Le[data_Le$delt_Le == "large_EcM",]

# Set factor order
data_Le$delt_Le <- factor(data_Le$delt_Le, levels = c("small_EcM", "large_EcM"))

# Merge with plot count
data_Le2 = data_Le %>% left_join(eco_plot_num, by = "ECOSUBCD")

# Create the recovery scatterplot
p2_1 = ggplot(data_Le2, aes(x = EcM_Le, y = AM_Le)) +
  geom_errorbar(aes(ymin = AM_Le - AM_Le_sd, ymax = AM_Le + AM_Le_sd),
                width = 0.06, color = "grey50", size = 0.4, alpha = 0.5) +
  geom_errorbarh(aes(xmin = EcM_Le - EcM_Le_sd, xmax = EcM_Le + EcM_Le_sd),
                 height = 0.06, color = "grey50", size = 0.4, alpha = 0.5) +
  geom_point(aes(fill = delt_Le, size = ECO_count),
             color = "black", shape = 21, alpha = 0.85) +
  scale_fill_manual(values = barCOLS, guide = "none") +
  scale_size_continuous(range = c(3, 12),
                        name = "Number of plots\nin ecoregion",
                        breaks = c(50, 1000, 2000, 3000),
                        guide = guide_legend(direction = "horizontal")) +
  scale_y_continuous(name = "Recovery of\nAM monospecific forests",
                     limits = c(1.5, 8.5), breaks = seq(2, 8, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(name = "Recovery of\nEcM monospecific forests",
                     limits = c(1.5, 8.5), breaks = seq(2, 8, 1),
                     expand = c(0, 0)) +
  geom_abline(slope = 1, linetype = "dashed", color = "black", size = 0.6, alpha = 0.8) +
  geom_polygon(aes(x = x, y = y), data = trsup_2, fill = "dodgerblue", alpha = 0.1) +
  geom_polygon(aes(x = x, y = y), data = trinf_2, fill = "mediumseagreen", alpha = 0.1) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.3), "cm"),
        legend.position = "bottom",
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 15, hjust = 0.5),
        legend.text = element_text(size = 15),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width  = unit(1.5, 'cm'),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.6),
        axis.text = element_text(size = 15),
        axis.text.y = element_text(hjust = 0.5, angle = 90),
        axis.ticks.length = unit(1, "mm"),
        strip.text = element_text(size = 12, family = "Arial", color = "black")) +
  annotate(geom = "text", x = 5, y = 8.2, label = Le_down, size = 6, color = "dodgerblue") +
  annotate(geom = "text", x = 5, y = 1.8, label = Le_up, size = 6, color = "mediumseagreen")

p2_1  # Print recovery plot
