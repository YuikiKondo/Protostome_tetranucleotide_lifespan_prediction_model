setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/01_lifespan_prediction_and_error/01_02_Fig2_weight")

# ==============================
# Positive/Negative WEIGHT SUMS by Region (all features)
# ==============================

sel_tag <- "ALL"   # label for output files

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Load
wdf <- read_csv("combined_model_weights_with_bio.csv", show_col_types = FALSE)

# Canonical region order -> pretty labels (match your other plots)
region_levels <- c("genome", "upstream", "exon", "intron", "downstream")
region_labels <- c("Genome", "Upstream", "Exon", "Intron", "Downstream")

# Filter and compute sums (no Selection_Frequency filter now)
region_weight_sums <- wdf %>%
  filter(!is.na(Average_Weight),
         Average_Weight != 0) %>%
  group_by(Region) %>%
  summarize(
    Positive_Weight_Sum = sum(Average_Weight[Average_Weight > 0], na.rm = TRUE),
    Negative_Weight_Sum = sum(Average_Weight[Average_Weight < 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(Region = region_levels,
           fill = list(Positive_Weight_Sum = 0, Negative_Weight_Sum = 0)) %>%
  mutate(
    Region = factor(Region, levels = region_levels, labels = region_labels)
  )

# Save the table
write_csv(region_weight_sums,
          paste0("weight_sums_by_region_", sel_tag, ".csv"))

# Long format for plotting
region_weight_sums_long <- region_weight_sums %>%
  pivot_longer(
    cols = c(Positive_Weight_Sum, Negative_Weight_Sum),
    names_to = "Polarity",
    values_to = "Weight_Sum"
  ) %>%
  mutate(
    Polarity = ifelse(Polarity == "Positive_Weight_Sum", "Positive", "Negative")
  )

# Plot
p_weight <- ggplot(region_weight_sums_long,
                   aes(x = Region, y = Weight_Sum, fill = Polarity)) +
  geom_col(position = "dodge", alpha = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey40") +
  scale_fill_manual(values = c("Positive" = "#d73027", "Negative" = "#4575b4")) +
  theme_minimal(base_size = 18) +
  labs(
    x = "Region",
    y = "Total weight sum",
    fill = "Polarity"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    axis.title.y = element_text(margin = margin(r = 15))
  )

ggsave(paste0("barplot_total_weight_sum_by_region_", sel_tag, ".png"),
       p_weight, width = 7, height = 5, dpi = 300)
print(p_weight)
