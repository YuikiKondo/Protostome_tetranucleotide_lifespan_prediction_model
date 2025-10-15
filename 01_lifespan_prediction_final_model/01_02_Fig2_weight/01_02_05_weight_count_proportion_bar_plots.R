
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# Tag to keep filenames consistent with the Selection_Frequency filter (> 0.1)
sel_tag <- "sel0.1"

# Helper: nice facet titles with preserved acronyms
pretty_trait <- function(x) {
  dplyr::recode(x,
    "AT_rich"        = "AT-rich",
    "CpG_including"  = "CpG-including",
    "CAA_including"  = "CAA-including",
    "ATGG_related"   = "ATGG-related",
    .default = x
  )
}


# ==============================
# 1) Total positive/negative COUNTS by Region (unchanged)
# ==============================

region_totals <- read_csv(paste0("count_region_totals_ALL_", sel_tag, ".csv"), show_col_types = FALSE)

region_totals_long <- region_totals %>%
  select(Region, Positive_Total_Region, Negative_Total_Region) %>%
  pivot_longer(cols = c(Positive_Total_Region, Negative_Total_Region),
               names_to = "Polarity", values_to = "Feature_Count") %>%
  mutate(
    Polarity = ifelse(Polarity == "Positive_Total_Region", "Positive", "Negative"),
    Region = factor(
      Region,
      levels = c("genome", "upstream", "exon", "intron", "downstream"),
      labels = c("Genome", "Upstream", "Exon", "Intron", "Downstream")
    )
  )

p_totals <- ggplot(region_totals_long, aes(x = Region, y = Feature_Count, fill = Polarity)) +
  geom_col(position = "dodge", alpha = 0.9) +
  scale_fill_manual(values = c("Positive" = "#d73027", "Negative" = "#4575b4")) +
  theme_minimal(base_size = 18) +
  labs(
    x = "Region",
    y = "Total feature count\n(sel.freq. > 0.1)",
    fill = "Polarity"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    axis.title.y = element_text(margin = margin(r = 15))
  )

ggsave(paste0("barplot_total_feature_counts_by_region_", sel_tag, ".png"),
       p_totals, width = 7, height = 5, dpi = 300)
print(p_totals)

# ==============================
# 2) Bar plots of COUNT proportions (Selection_Frequency > 0.1)
#    — show only selected traits, add Possible_Tetranucleotide_Number in titles,
#      and Expected_Prop as dotted line
# ==============================

props <- read_csv(paste0("count_proportions_by_trait_and_region_", sel_tag, ".csv"),
                  show_col_types = FALSE)

selected_traits <- c("ATGG_related", "CpG_including", "CAA_including", "AT_rich")

# Build a lookup for possible tetranucleotide numbers and expected proportions
trait_info <- props %>%
  distinct(Trait, Possible_Tetranucleotide_Number, Expected_Prop) %>%
  filter(Trait %in% selected_traits)

props_long <- props %>%
  filter(Trait %in% selected_traits) %>%
  mutate(Trait = factor(Trait, levels = selected_traits)) %>%
  select(Trait, Region, Positive_Prop, Negative_Prop) %>%
  pivot_longer(cols = c(Positive_Prop, Negative_Prop),
               names_to = "Polarity", values_to = "Proportion") %>%
  mutate(
    Polarity = ifelse(Polarity == "Positive_Prop", "Positive", "Negative"),
    Region = factor(
      Region,
      levels = c("genome", "upstream", "exon", "intron", "downstream"),
      labels = c("Genome", "Upstream", "Exon", "Intron", "Downstream")
    )
  ) %>%
  left_join(trait_info, by = "Trait") %>%
  mutate(
    Trait_pretty = paste0(
      pretty_trait(as.character(Trait)),
      " (", Possible_Tetranucleotide_Number, ")"
    )
  )

p <- ggplot(props_long, aes(x = Region, y = Proportion, fill = Polarity)) +
  geom_col(position = "dodge", alpha = 0.9) +
  facet_wrap(~ Trait_pretty, nrow = 1, scales = "fixed") +
  scale_fill_manual(values = c("Positive" = "#d73027", "Negative" = "#4575b4")) +
  scale_y_continuous(limits = c(0, 1)) +
  # dotted line for expected proportion
  geom_hline(aes(yintercept = Expected_Prop),
             data = distinct(props_long, Trait_pretty, Expected_Prop),
             linetype = "dotted", color = "black", size = 1) +
  theme_minimal(base_size = 18) +
  labs(
    x = "Region",
    y = "Proportion of features\n(sel.freq. > 0.1)",
    fill = "Polarity"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold"),
    legend.position = "right",
    axis.title.y = element_text(margin = margin(r = 15))
  )

ggsave(paste0("barplot_count_proportions_by_trait_and_region_", sel_tag, "_4traits_with_numbers_expected.png"),
       p, width = 14, height = 5, dpi = 300)
print(p)
