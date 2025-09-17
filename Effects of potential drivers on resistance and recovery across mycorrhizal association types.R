############################################################
## Dependencies (load before running):
##   library(dplyr); library(tidyr); library(tibble); library(ggplot2)
##   library(nlme); library(car); library(forcats)
##   (Optional) library(performance) for multicollinearity diagnostics on mixed models
############################################################

### Load the dataset Arbuscularmycorrhizalforests, Ectomycorrhizalforests, Mixedmycorrhizalforests ###

###### Mixedness index ########################################################
# Define a function that assigns the highest “mixedness” when EcM = 50% and
# linearly decreases toward 0 as EcM approaches 0% or 100%.
# - x: EcM percentage in [0, 100]
# - max_value: the peak value at exactly 50% EcM (default 100)
# Output: a number in [0, max_value], triangular around x = 50.
calculate_mixed <- function(x, max_value = 100) {
  mixed_value <- max_value * (1 - abs(x - 50) / 50)
  return(mixed_value)
}

# Apply the mixedness function to each dataset’s EcM percentage.
# NOTE: sapply() is fine; x is numeric. If EcM_percentage can contain NA,
# calculate_mixed handles it implicitly as arithmetic with NA returns NA.
data_all$mixed_percentage   <- sapply(data_all$EcM_percentage,   calculate_mixed)
data_AM$mixed_percentage    <- sapply(data_AM$EcM_percentage,    calculate_mixed)
data_EcM$mixed_percentage   <- sapply(data_EcM$EcM_percentage,   calculate_mixed)
data_EcM_AM$mixed_percentage<- sapply(data_EcM_AM$EcM_percentage,calculate_mixed)

###### Standardization (z-scoring) ###########################################
# Standardize selected columns (mean 0, SD 1) to make coefficients comparable.
# NOTE: base::scale() returns a matrix with attributes; assigning it to a
# data.frame column is fine. If you need a pure numeric vector, use as.numeric(scale(...)).
# Ensure the listed columns exist and are numeric.
data_all$species_richness_log <- scale(data_all$species_richness_log)
data_all$mixed_percentage     <- scale(data_all$mixed_percentage)
data_all$EcM_percentage       <- scale(data_all$EcM_percentage)
data_all$STDAGE               <- scale(data_all$STDAGE)
data_all$SLOPE                <- scale(data_all$SLOPE)
data_all$ELEV                 <- scale(data_all$ELEV)

data_all$ph                   <- scale(data_all$ph)
data_all$totaln               <- scale(data_all$totaln)
data_all$bulk                 <- scale(data_all$bulk)
data_all$clay                 <- scale(data_all$clay)

data_all$pre                  <- scale(data_all$pre)
data_all$tem                  <- scale(data_all$tem)

# Repeat standardization for AM-only dataset
data_AM$species_richness_log <- scale(data_AM$species_richness_log)
data_AM$mixed_percentage     <- scale(data_AM$mixed_percentage)
data_AM$EcM_percentage       <- scale(data_AM$EcM_percentage)
data_AM$STDAGE               <- scale(data_AM$STDAGE)
data_AM$SLOPE                <- scale(data_AM$SLOPE)
data_AM$ELEV                 <- scale(data_AM$ELEV)

data_AM$ph                   <- scale(data_AM$ph)
data_AM$totaln               <- scale(data_AM$totaln)
data_AM$bulk                 <- scale(data_AM$bulk)
data_AM$clay                 <- scale(data_AM$clay)

data_AM$pre                  <- scale(data_AM$pre)
data_AM$tem                  <- scale(data_AM$tem)

# Repeat for EcM-only dataset
data_EcM$species_richness_log <- scale(data_EcM$species_richness_log)
data_EcM$mixed_percentage     <- scale(data_EcM$mixed_percentage)
data_EcM$EcM_percentage       <- scale(data_EcM$EcM_percentage)
data_EcM$STDAGE               <- scale(data_EcM$STDAGE)
data_EcM$SLOPE                <- scale(data_EcM$SLOPE)
data_EcM$ELEV                 <- scale(data_EcM$ELEV)

data_EcM$ph                   <- scale(data_EcM$ph)
data_EcM$totaln               <- scale(data_EcM$totaln)
data_EcM$bulk                 <- scale(data_EcM$bulk)
data_EcM$clay                 <- scale(data_EcM$clay)

data_EcM$pre                  <- scale(data_EcM$pre)
data_EcM$tem                  <- scale(data_EcM$tem)

# Repeat for mixed EcM–AM dataset
data_EcM_AM$species_richness_log <- scale(data_EcM_AM$species_richness_log)
data_EcM_AM$mixed_percentage     <- scale(data_EcM_AM$mixed_percentage)
data_EcM_AM$EcM_percentage       <- scale(data_EcM_AM$EcM_percentage)
data_EcM_AM$STDAGE               <- scale(data_EcM_AM$STDAGE)
data_EcM_AM$SLOPE                <- scale(data_EcM_AM$SLOPE)
data_EcM_AM$ELEV                 <- scale(data_EcM_AM$ELEV)

data_EcM_AM$ph                   <- scale(data_EcM_AM$ph)
data_EcM_AM$totaln               <- scale(data_EcM_AM$totaln)
data_EcM_AM$bulk                 <- scale(data_EcM_AM$bulk)
data_EcM_AM$clay                 <- scale(data_EcM_AM$clay)

data_EcM_AM$pre                  <- scale(data_EcM_AM$pre)
data_EcM_AM$tem                  <- scale(data_EcM_AM$tem)

##############################################################################
## Resistance (Rs) models
##############################################################################
library(car)   # for vif()
# NOTE: Models below use nlme::lme; ensure library(nlme) is attached beforehand.
# Random intercept: ECOSUBCD (ecoregion/subsection)
# Spatial correlation: exponential correlation (corExp) on geographic coordinates LON, LAT
# Response: Rs_log; Predictors standardized above.

# ----- Total dataset model (includes mixed and EcM percentages) -----
model_data_Total <- lme(
  Rs_log ~ species_richness_log + mixed_percentage + EcM_percentage +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,                       # random intercept per ECOSUBCD
  correlation = corExp(form = ~ LON + LAT),      # exponential spatial correlation
  data = data_all
)
dim(data_all); vif(model_data_Total)  # NOTE: car::vif() is designed for lm/glm; for lme this may not work in some R versions.

# Tidy coefficient table: drop intercept, strip "scale()" from names, add CI = beta ± 1.96*SE
model_Total <- summary(model_data_Total)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "Total") %>%
  mutate(conf.low = Value - 1.96 * Std.Error,
         conf.high = Value + 1.96 * Std.Error)

# ----- AM-only model (no EcM/mixed % terms by design) -----
model_data_AM <- lme(
  Rs_log ~ species_richness_log +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,
  correlation = corExp(form = ~ LON + LAT),
  data = filtered_data_AM
)
dim(filtered_data_AM); vif(model_data_AM)

model_AM <- summary(model_data_AM)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "AM") %>%
  mutate(conf.low = Value - 1.96 * Std.Error,
         conf.high = Value + 1.96 * Std.Error)

# ----- EcM-only model -----
model_data_EcM <- lme(
  Rs_log ~ species_richness_log +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,
  correlation = corExp(form = ~ LON + LAT),
  data = filtered_data_EcM
)
dim(filtered_data_EcM); vif(model_data_EcM)

model_EcM <- summary(model_data_EcM)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "EcM") %>%
  mutate(conf.low = Value - 1.96 * Std.Error,
         conf.high = Value + 1.96 * Std.Error)

# ----- Mixed EcM–AM model -----
model_data_EcM_AM <- lme(
  Rs_log ~ species_richness_log +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,
  correlation = corExp(form = ~ LON + LAT),
  data = filtered_data_EcM_AM
)
dim(filtered_data_EcM_AM); vif(model_data_EcM_AM)

model_EcM_AM <- summary(model_data_EcM_AM)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "Mixed EcM and AM") %>%
  mutate(conf.low = Value - 1.96 * Std.Error,
         conf.high = Value + 1.96 * Std.Error)

# Combine all Rs coefficient tables and classify predictors into driver groups.
model_total_Rs <- rbind(model_Total, model_EcM, model_AM, model_EcM_AM) %>%
  mutate(
    group = case_when(
      factor %in% c("LON", "LAT", "ELEV", "SLOPE") ~ "Geographical attributes",
      factor %in% c("mixed_percentage", "EcM_percentage") ~ "Mycorrhizas",
      factor %in% c("ph", "totaln", "bulk", "clay") ~ "Soil properties",
      factor %in% c("pre", "tem") ~ "Climatic conditions",
      factor %in% c("STDAGE", "tree_count") ~ "Stand age",
      factor %in% c("species_richness_log") ~ "Species richness",
      TRUE ~ "Other"
    ),
    # Human-readable factor labels for plotting
    factor = case_when(
      factor == "species_richness_log" ~ "Species richness",
      factor == "STDAGE" ~ "Stand age",
      factor == "tree_count" ~ "Tree count",
      factor == "mixed_percentage" ~ "Mixed mycorrhizal associations",
      factor == "EcM_percentage" ~ "EcM percentage",
      factor == "ph" ~ "pH",
      factor == "totaln" ~ "Total nitrogen content",
      factor == "bulk" ~ "Bulk density",
      factor == "clay" ~ "Clay content",
      factor == "pre" ~ "Mean annual precipitation",
      factor == "tem" ~ "Mean annual temperature",
      factor == "LON" ~ "Longitude",
      factor == "LAT" ~ "Latitude",
      factor == "SLOPE" ~ "Slope",
      factor == "ELEV" ~ "Elevation",
      TRUE ~ factor
    )
  )

# Control the display order of groups, factors, and mycorrhizal types in plots.
model_total_Rs$group <- factor(
  model_total_Rs$group,
  levels = c("Climatic conditions","Soil properties","Geographical attributes",
             "Stand age","Species richness","Mycorrhizas")
)
model_total_Rs$factor <- factor(
  model_total_Rs$factor,
  levels = c("EcM percentage","Mixed mycorrhizal associations","Species richness","Stand age",
             "Mean annual precipitation","Mean annual temperature",
             "pH","Bulk density","Total nitrogen content","Clay content",
             "Elevation","Slope")
)
model_total_Rs$Mycorrhizal.type <- factor(
  model_total_Rs$Mycorrhizal.type,
  levels = c("Total","EcM","AM","Mixed EcM and AM")
)

library(forcats)

# Colors for dots/bars (consistent across driver groups)
dotCOLS <- c("#99C3AF","#C9E8C9","#DFEBD6",
             "#9CB9D7","#C6ACD9","#F4B7AD")
barCOLS <- c("#006837","#78C679","#aecc99",
             "#08519C","#7030A0","#E34A33")

# Plot p1: standardized coefficients for Rs (point + 95% CI error bars), flipped coordinates.
p1 <- ggplot(
  model_total_Rs,
  aes(x = factor, y = Value, ymin = conf.low, ymax = conf.high,
      fill = group, color = group, shape = Mycorrhizal.type)
) +
  geom_hline(yintercept = 0, lty = 2) +  # reference line at zero effect
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.8)) +
  geom_point(size = 3, position = position_dodge(width = 0.8)) +
  scale_fill_manual(name = "Potential driver type", values = barCOLS) +
  scale_color_manual(name = "Potential driver type", values = dotCOLS) +
  scale_shape_manual(name = "Mycorrhizal type", values = c(22, 21, 23, 24)) +
  labs(title = expression(Resistance ~ (italic(R)[s]))) +
  xlab(NULL) +
  theme_bw() +
  scale_y_continuous(
    name   = "Standardized regression coefficients",
    limits = c(-0.53, 0.28),
    breaks = seq(-0.5, 0.25, 0.25),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_flip() +  # horizontal layout: factors on y-axis
  theme(
    plot.margin       = unit(c(0.5,0.2,0.1,0.2), "cm"),
    axis.title.x      = element_text(size = 16),
    axis.title.y      = element_text(size = 16),
    legend.position   = "none",
    plot.title        = element_text(hjust = 0.5, size = 16),
    legend.title      = element_text(size = 15),
    legend.text       = element_text(size = 15),
    legend.key.height = unit(0.25, 'cm'),
    legend.key.width  = unit(0.75, 'cm'),
    axis.text.x       = element_text(size = 15),
    axis.text.y       = element_text(size = 15, color = "black"),
    axis.ticks.length = unit(1, "mm"),
    strip.text        = element_text(size = 16, family = "Arial", color = "black"),
    legend.spacing.x  = unit(2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_rect(color = "grey40", fill = "white")
  ) +
  guides(
    fill  = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
    color = guide_legend(title.position = "top", title.hjust = 0.5,
                         override.aes = list(fill = NA, shape = NA), nrow = 2),
    shape = guide_legend(title.position = "top", title.hjust = 0.5,
                         override.aes = list(fill = "grey50"), nrow = 2)
  )
p1

# Aggregate absolute standardized coefficients by driver group per mycorrhizal type.
# Interpreted as group-level “relative importance” (sum of |beta|).
model_total_Rs_group <- rbind(model_Total, model_EcM, model_AM, model_EcM_AM) %>%
  mutate(
    group = case_when(
      factor %in% c("LON", "LAT", "ELEV", "SLOPE") ~ "Geographical attributes",
      factor %in% c("mixed_percentage", "EcM_percentage") ~ "Mycorrhizas",
      factor %in% c("ph", "totaln", "bulk", "clay") ~ "Soil properties",
      factor %in% c("pre", "tem") ~ "Climatic conditions",
      factor %in% c("STDAGE", "tree_count") ~ "Stand age",
      factor %in% c("species_richness_log") ~ "Species richness",
      TRUE ~ "Other"
    ),
    factor = case_when(
      factor == "species_richness_log" ~ "Species richness",
      factor == "STDAGE" ~ "Stand age",
      factor == "tree_count" ~ "Tree count",
      factor == "mixed_percentage" ~ "Mixed mycorrhizal associations",
      factor == "EcM_percentage" ~ "EcM percentage",
      factor == "ph" ~ "pH",
      factor == "totaln" ~ "Total nitrogen content",
      factor == "bulk" ~ "Bulk density",
      factor == "clay" ~ "Clay content",
      factor == "pre" ~ "Mean annual precipitation",
      factor == "tem" ~ "Mean annual temperature",
      factor == "LON" ~ "Longitude",
      factor == "LAT" ~ "Latitude",
      factor == "SLOPE" ~ "Slope",
      factor == "ELEV" ~ "Elevation",
      TRUE ~ factor
    )
  ) %>%
  group_by(group, Mycorrhizal.type) %>%
  summarise(total_value = sum(abs(Value)))
print(model_total_Rs_group)
range(model_total_Rs_group$total_value)

model_total_Rs_group$group <- factor(
  model_total_Rs_group$group,
  levels = c("Climatic conditions","Soil properties","Geographical attributes",
             "Stand age","Species richness","Mycorrhizas")
)

# Bar plot p2: group-level relative importance for Rs by mycorrhizal type (facets)
p2 <- ggplot(model_total_Rs_group, aes(x = group, y = total_value, fill = group, color = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", alpha = 0.95) +
  facet_wrap(~ fct_relevel(Mycorrhizal.type, "Total", "EcM", "AM", "Mixed EcM and AM"),
             ncol = 1, strip.position = "right") +
  scale_y_continuous(breaks = seq(0, 0.6, 0.3), limits = c(0, 0.6),
                     expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = NULL, y = expression(The~relative~importance~of~resistance~(italic(R)[s]))) +
  scale_fill_manual(values = c(
    "Climatic conditions" = "#99C3AF",
    "Soil properties"     = "#C9E8C9",
    "Geographical attributes" = "#DFEBD6",
    "Stand age"           = "#9CB9D7",
    "Species richness"    = "#C6ACD9",
    "Mycorrhizas"         = "#F4B7AD"
  )) +
  scale_color_manual(values = c(
    "Climatic conditions" = "#99C3AF",
    "Soil properties"     = "#C9E8C9",
    "Geographical attributes" = "#DFEBD6",
    "Stand age"           = "#9CB9D7",
    "Species richness"    = "#C6ACD9",
    "Mycorrhizas"         = "#F4B7AD"
  )) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 11),
    strip.background.x = element_rect(color = "#EEEEEF", fill = "#EEEEEF"),
    strip.background.y = element_rect(color = "#EEEEEF", fill = "#EEEEEF"),
    strip.placement = "outside",
    panel.spacing.x = unit(0.3, "cm"),
    panel.spacing.y = unit(0.3, "cm"),
    plot.margin     = unit(c(1.5,0.2,0.1,0.5), "cm"),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    axis.text.x     = element_text(size = 14, vjust = 1, hjust = 1, angle = 45, color = "black"),
    axis.text.y     = element_text(hjust = 0.5, size = 14, angle = 90),
    plot.title      = element_text(hjust = 0.5, size = 18),
    legend.position = "none",
    legend.title    = element_text(hjust = 0.5, size = 13.5, family = "Arial", color = "black"),
    legend.text     = element_text(size = 13.5, family = "Arial", color = "black"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width  = unit(1.25, 'cm'),
    legend.key         = element_rect(fill = "white"),
    legend.background  = element_blank()
  )
p2

##############################################################################
## Recovery (Le) models
##############################################################################
library(car)

# ----- Total dataset model for Le_log (includes mixed and EcM percentages) -----
model_data_Total_Le <- lme(
  Le_log ~ species_richness_log + mixed_percentage + EcM_percentage +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,
  correlation = corExp(form = ~ LON + LAT),
  data = data_all
)
dim(data_all)

model_Total_Le <- summary(model_data_Total_Le)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "Total") %>%
  mutate(conf.low = Value - 1.96 * Std.Error,
         conf.high = Value + 1.96 * Std.Error)

# ----- AM-only Le model (subset filtered by ECOSUBCD with >30 plots if wanted) -----
# The block below identifies ECOSUBCD with at least 31 distinct plot_ID (n > 30)
# and prepares a filtered dataset; the model just below still uses filtered_data_AM.
valid_groups <- filtered_data_AM %>%
  group_by(ECOSUBCD) %>%
  summarise(n = n_distinct(plot_ID)) %>%
  filter(n > 30)

filtered_data_AM_filtered <- filtered_data_AM %>%
  semi_join(valid_groups, by = "ECOSUBCD")

model_data_AM_Le <- lme(
  Le_log ~ species_richness_log +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,
  correlation = corExp(form = ~ LON + LAT),
  data = filtered_data_AM
)
summary(model_data_AM_Le)
dim(filtered_data_AM)

model_AM_Le <- summary(model_data_AM_Le)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "AM") %>%
  # NOTE: Here you use ±1*SE (approx 68% CI), not 1.96*SE. Kept as-is intentionally.
  mutate(conf.low = Value - 1 * Std.Error,
         conf.high = Value + 1 * Std.Error)

# ----- EcM-only Le model -----
model_data_EcM_Le <- lme(
  Le_log ~ species_richness_log +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,
  correlation = corExp(form = ~ LON + LAT),
  data = filtered_data_EcM
)
summary(model_data_EcM_Le)
vif(model_data_EcM)    # NOTE: calling VIF on the Rs EcM model; left unchanged.
dim(data_EcM)

model_EcM_Le <- summary(model_data_EcM_Le)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "EcM") %>%
  mutate(conf.low = Value - 1.96 * Std.Error,
         conf.high = Value + 1.96 * Std.Error)

# ----- Mixed EcM–AM Le model (no explicit spatial correlation here; kept as provided) -----
model_data_EcM_AM_Le <- lme(
  Le_log ~ species_richness_log +
    STDAGE + SLOPE + ELEV +
    ph + totaln + bulk + clay +
    pre + tem,
  random = ~ 1 | ECOSUBCD,
  data = filtered_data_EcM_AM
)
summary(model_data_EcM_AM_Le)
dim(data_EcM_AM)

model_EcM_AM_Le <- summary(model_data_EcM_AM_Le)$tTable %>%
  as.data.frame() %>%
  rownames_to_column(var = "factor") %>%
  filter(factor != "(Intercept)") %>%
  mutate(factor = gsub("scale\\(|\\)", "", factor)) %>%
  mutate(Mycorrhizal.type = "Mixed EcM and AM") %>%
  mutate(conf.low = Value - 1.96 * Std.Error,
         conf.high = Value + 1.96 * Std.Error)

# Combine Le coefficients and classify into groups (includes “Resistance” as a predictor where present).
model_total_Le <- rbind(model_Total_Le, model_EcM_Le, model_AM_Le, model_EcM_AM_Le) %>%
  mutate(
    group = case_when(
      factor %in% c("LON", "LAT", "ELEV", "SLOPE") ~ "Geographical attributes",
      factor %in% c("mixed_percentage", "EcM_percentage") ~ "Mycorrhizas",
      factor %in% c("ph", "totaln", "bulk", "clay") ~ "Soil properties",
      factor %in% c("pre", "tem") ~ "Climatic conditions",
      factor %in% c("STDAGE", "tree_count") ~ "Stand structure",
      factor %in% c("species_richness_log") ~ "Species richness",
      factor %in% c("Rs_log") ~ "Resistance",
      TRUE ~ "Other"
    ),
    factor = case_when(
      factor == "species_richness_log" ~ "Species richness",
      factor == "STDAGE" ~ "Stand age",
      factor == "Rs_log" ~ "Resistance",
      factor == "mixed_percentage" ~ "Mixed mycorrhizal associations",
      factor == "EcM_percentage" ~ "EcM percentage",
      factor == "ph" ~ "pH",
      factor == "totaln" ~ "Total nitrogen content",
      factor == "bulk" ~ "Bulk density",
      factor == "clay" ~ "Clay content",
      factor == "pre" ~ "Mean annual precipitation",
      factor == "tem" ~ "Mean annual temperature",
      factor == "LON" ~ "Longitude",
      factor == "LAT" ~ "Latitude",
      factor == "SLOPE" ~ "Slope",
      factor == "ELEV" ~ "Elevation",
      TRUE ~ factor
    )
  )

# Order for plotting
model_total_Le$group <- factor(
  model_total_Le$group,
  levels = c("Climatic conditions","Soil properties","Geographical attributes",
             "Stand structure","Species richness","Mycorrhizas","Resistance")
)
model_total_Le$factor <- factor(
  model_total_Le$factor,
  levels = c("Resistance","EcM percentage","Mixed mycorrhizal associations",
             "Species richness","Stand age",
             "Mean annual precipitation","Mean annual temperature",
             "pH","Bulk density","Total nitrogen content","Clay content",
             "Elevation","Slope")
)
model_total_Le$Mycorrhizal.type <- factor(
  model_total_Le$Mycorrhizal.type,
  levels = c("Total","EcM","AM","Mixed EcM and AM")
)

library(forcats)

# Distinct color palettes for Le plot (includes “Resistance” group)
dotCOLS_Le <- c("#99C3AF","#C9E8C9","#DFEBD6",
                "#9CB9D7","#C6ACD9","#F4B7AD","#D1D17F")
barCOLS_Le <- c("#006837","#78C679","#aecc99",
                "#08519C","#7030A0","#E34A33","#a4a400")

# Quick range checks for CI bounds (useful to set y-limits)
max(model_total_Le$conf.high)
min(model_total_Le$conf.low)

# Plot p3: standardized coefficients for Le (point + CI), flipped coordinates.
p3 <- ggplot(
  model_total_Le,
  aes(x = factor, y = Value, ymin = conf.low, ymax = conf.high,
      fill = group, color = group, shape = Mycorrhizal.type)
) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.2, position = position_dodge(width = 0.8)) +
  geom_point(size = 3, position = position_dodge(width = 0.8)) +
  scale_fill_manual(name = "Potential driver type", values = barCOLS_Le) +
  scale_color_manual(name = "Potential driver type", values = dotCOLS_Le) +
  scale_shape_manual(name = "Mycorrhizal type", values = c(22, 21, 23, 24)) +
  labs(title = expression(Recovery ~ (italic(R)[c]))) +
  xlab(NULL) +
  theme_bw() +
  scale_y_continuous(
    name   = "Standardized regression coefficients",
    limits = c(-0.32, 0.32),
    breaks = seq(-0.30, 0.30, 0.15),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_flip() +
  theme(
    plot.margin       = unit(c(0.5,0.2,0.1,0.2), "cm"),
    axis.title.x      = element_text(size = 16),
    axis.title.y      = element_text(size = 16),
    legend.position   = "bottom",
    plot.title        = element_text(hjust = 0.5, size = 16),
    legend.title      = element_text(size = 15),
    legend.text       = element_text(size = 15),
    legend.key.height = unit(0.25, 'cm'),
    legend.key.width  = unit(0.75, 'cm'),
    axis.text.x       = element_text(size = 15),
    axis.text.y       = element_text(size = 15, color = "black"),
    axis.ticks.length = unit(1, "mm"),
    strip.text        = element_text(size = 16, family = "Arial", color = "black"),
    legend.spacing.x  = unit(2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_rect(color = "grey40", fill = "white")
  ) +
  guides(
    fill  = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2),
    color = guide_legend(title.position = "top", title.hjust = 0.5,
                         override.aes = list(fill = NA, shape = NA), nrow = 2),
    shape = guide_legend(title.position = "top", title.hjust = 0.5,
                         override.aes = list(fill = "grey50"), nrow = 2)
  )
p3

# Aggregate absolute standardized coefficients for Le by group and mycorrhizal type.
model_total_Le_group <- rbind(model_Total_Le, model_EcM_Le, model_AM_Le, model_EcM_AM_Le) %>%
  mutate(
    group = case_when(
      factor %in% c("LON", "LAT", "ELEV", "SLOPE") ~ "Geographical attributes",
      factor %in% c("mixed_percentage", "EcM_percentage") ~ "Mycorrhizas",
      factor %in% c("ph", "totaln", "bulk", "clay") ~ "Soil properties",
      factor %in% c("pre", "tem") ~ "Climatic conditions",
      factor %in% c("STDAGE", "tree_count") ~ "Stand age",
      factor %in% c("species_richness_log") ~ "Species richness",
      factor %in% c("Rs_log") ~ "Resistance",
      TRUE ~ "Other"
    ),
    factor = case_when(
      factor == "species_richness_log" ~ "Species richness",
      factor == "STDAGE" ~ "Stand age",
      factor == "Rs_log" ~ "Resistance",
      factor == "mixed_percentage" ~ "Mixed mycorrhizal associations",
      factor == "EcM_percentage" ~ "EcM percentage",
      factor == "ph" ~ "pH",
      factor == "totaln" ~ "Total nitrogen content",
      factor == "bulk" ~ "Bulk density",
      factor == "clay" ~ "Clay content",
      factor == "pre" ~ "Mean annual precipitation",
      factor == "tem" ~ "Mean annual temperature",
      factor == "LON" ~ "Longitude",
      factor == "LAT" ~ "Latitude",
      factor == "SLOPE" ~ "Slope",
      factor == "ELEV" ~ "Elevation",
      TRUE ~ factor
    )
  ) %>%
  group_by(group, Mycorrhizal.type) %>%
  summarise(total_value = sum(abs(Value)))
print(model_total_Le_group)
range(model_total_Le_group$total_value)

model_total_Le_group$group <- factor(
  model_total_Le_group$group,
  levels = c("Climatic conditions","Soil properties","Geographical attributes",
             "Stand age","Species richness","Mycorrhizas","Resistance")
)

# Bar plot p4: group-level relative importance for Le
p4 <- ggplot(model_total_Le_group, aes(x = group, y = total_value, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", alpha = 0.95) +
  scale_y_continuous(breaks = seq(0, 0.6, 0.3), limits = c(0, 0.6),
                     expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap(~ fct_relevel(Mycorrhizal.type, "Total", "EcM", "AM", "Mixed EcM and AM"),
             ncol = 1, strip.position = "right") +
  labs(x = NULL, y = expression(The~relative~importance~of~recovery~(italic(R)[c]))) +
  scale_fill_manual(values = c(
    "Climatic conditions"     = "#99C3AF",
    "Soil properties"         = "#C9E8C9",
    "Geographical attributes" = "#DFEBD6",
    "Resistance"              = "#D1D17F",
    "Stand age"               = "#9CB9D7",
    "Species richness"        = "#C6ACD9",
    "Mycorrhizas"             = "#F4B7AD"
  )) +
  theme_bw() +
  theme(
    strip.text.x    = element_text(size = 15),
    strip.text.y    = element_text(size = 11),
    strip.background.x = element_rect(color = "#EEEEEF", fill = "#EEEEEF"),
    strip.background.y = element_rect(color = "#EEEEEF", fill = "#EEEEEF"),
    strip.placement = "outside",
    panel.spacing.x = unit(0.3, "cm"),
    panel.spacing.y = unit(0.3, "cm"),
    plot.margin     = unit(c(1.5,0.2,0.1,0.5), "cm"),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    axis.text.x     = element_text(size = 14, vjust = 1, hjust = 1, angle = 45, color = "black"),
    axis.text.y     = element_text(hjust = 0.5, size = 14, angle = 90),
    plot.title      = element_text(hjust = 0.5, size = 18),
    legend.position = "none",
    legend.title    = element_text(hjust = 0.5, size = 13.5, family = "Arial", color = "black"),
    legend.text     = element_text(size = 13.5, family = "Arial", color = "black"),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width  = unit(1.25, 'cm'),
    legend.key         = element_rect(fill = "white"),
    legend.background  = element_blank()
  )