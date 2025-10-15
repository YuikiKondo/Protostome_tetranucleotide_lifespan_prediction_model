# =======================
# Setup
# =======================
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/05_mechanical_index_parameter_and_weight_correlation/mechanics_scatterplots_fdr_by_region_GLOBAL")

library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tools)
library(stats)
library(stringr)

# ---------------------------------
# Helpers
# ---------------------------------
clean_region <- function(x) {
  x <- tolower(x)
  x <- gsub("upstream2", "upstream", x)
  x <- gsub("downstream2", "downstream", x)
  x <- toTitleCase(x)
  x
}

# Find & rename a column if any of the aliases exist
rename_if_present <- function(df, aliases, new){
  hit <- intersect(names(df), aliases)
  if (length(hit) == 1) df <- dplyr::rename(df, !!new := !!sym(hit)) else if (length(hit) > 1) {
    # if multiple matches, prefer the first in aliases order
    keep <- hit[match(hit, aliases)]
    df <- dplyr::rename(df, !!new := !!sym(keep[1]))
  }
  df
}

# Make Index values comparable across files
normalize_index <- function(x){
  x %>%
    as.character() %>%
    sub("_z$", "", ., ignore.case = TRUE) # remove trailing "_z" if present
}

# =======================
# Load data
# =======================
spearman_df <- read_csv("spearman_p_and_GLOBAL_FDR_summary_ALL_REGIONS.csv", show_col_types = FALSE)
pearson_df  <- read_csv("pearson_correlation_results_DNA_mechanics_and_lifespan_global_FDR.csv", show_col_types = FALSE)

# -----------------------
# Harmonize columns
# -----------------------
# spearman_df is expected to have Region / Index / Spearman_rho (but make robust)
spearman_df <- spearman_df %>%
  rename_if_present(c("Region","region","Region_name"), "Region") %>%
  rename_if_present(c("Index","Mechanics_Index","mechanics_index","MechanicsIndex"), "Index") %>%
  rename_if_present(c("Spearman_rho","rho","SpearmanRho","spearman_rho"), "Spearman_rho") %>%
  mutate(
    Region = clean_region(Region),
    Index  = normalize_index(Index)
  )

# pearson_df may use different names; map flexibly
pearson_df <- pearson_df %>%
  rename_if_present(c("Region","region","Region_name"), "Region") %>%
  rename_if_present(c("Index","Mechanics_Index","mechanics_index","MechanicsIndex"), "Index") %>%
  # Pearson r column could be "Pearson_r", "r", "Correlation_r"
  rename_if_present(c("Pearson_r","r","Correlation_r","pearson_r"), "Pearson_r") %>%
  # (Optional meta columns that might exist; keep if you like)
  mutate(
    Region = clean_region(Region),
    Index  = normalize_index(Index)
  )

# Safety checks
stopifnot(all(c("Region","Index","Spearman_rho") %in% names(spearman_df)))
stopifnot(all(c("Region","Index","Pearson_r") %in% names(pearson_df)))

# =======================
# Merge on Region + Index
# =======================
merged <- inner_join(
  spearman_df %>% select(Region, Index, Spearman_rho),
  pearson_df  %>% select(Region, Index, Pearson_r),
  by = c("Region", "Index")
)

cat("Merged rows:", nrow(merged), "\n")
print(head(merged))

# =======================
# Correlation per region (Pearson)
# =======================
stats_df <- merged %>%
  group_by(Region) %>%
  summarise(
    test = list(cor.test(Pearson_r, Spearman_rho, method = "pearson")),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    r  = unname(test$estimate),
    p  = test$p.value
  ) %>%
  ungroup() %>%
  mutate(FDR_adjusted_p = p.adjust(p, method = "BH")) %>%
  select(Region, r, p, FDR_adjusted_p) %>%
  arrange(Region)

print(stats_df)

# =======================
# Annotation for facets
# =======================
annot_df <- stats_df %>%
  mutate(
    label = paste0(
      "r = ", sprintf("%.2f", r),
      "\nFDR p = ", format.pval(FDR_adjusted_p, digits = 2, eps = 1e-3)
    ),
    x = -0.55,   # left edge of x-axis
    y =  0.28,   # top of y-axis
    hjust = 0,
    vjust = 1
  )

# =======================
# Plot
# =======================
p <- ggplot(merged, aes(x = Pearson_r, y = Spearman_rho)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", color = "red") +
  facet_wrap(~ Region) +
  scale_x_continuous(limits = c(-0.55, 0.5)) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  labs(
    x = "Pearson's r (ln(lifespan days) ~ raw mechanical index)",
    y = "Spearman's rho (weights ~ mechanical index parameter)"
  ) +
  geom_text(
    data = annot_df,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(face = "bold"))

# Save PNG
ggsave(
  filename = "scatter_lifespan_vs_weights_by_region.png",
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

# =======================
# Save stats table
# =======================
write_csv(stats_df, "lifespan_vs_weight_scatter_stats_by_region.csv")
