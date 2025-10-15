# ====== Set working directory ======
setwd("/Users/yuikikondo/Desktop/Github_code/protostome_lifespan_prediction_model/03_variance_decomposition")

# ====== Config ======
csv_path <- "variance_partitioning_T_M_P_G_atoms.csv"

# Output names
out_atoms  <- "variance_partitioning_R2_atoms_sorted.png"
out_totals <- "variance_partitioning_R2_totals_by_set_bars.png"

# Colors
HILITE_BLUE <- "#1f78b4"
GREY_OTHER  <- "grey70"

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(ggplot2); library(tidyr)
})

# ====== Load ======
df <- read_csv(csv_path, show_col_types = FALSE)

# Expect columns: Atom, R2_amount, Percent_of_R2_TMPG
stopifnot(all(c("Atom","R2_amount","Percent_of_R2_TMPG") %in% names(df)))

# ====== Clip tiny numerical negatives ======
EPS <- 1e-9
df <- df |>
  mutate(
    R2_amount = ifelse(R2_amount < 0 & R2_amount > -EPS, 0, R2_amount),
    Percent_of_R2_TMPG = ifelse(Percent_of_R2_TMPG < 0 & Percent_of_R2_TMPG > -EPS, 0, Percent_of_R2_TMPG)
  )

# ------------------------------------------------------------------------------
# 1) Atom-sorted bar plot, highlight all T-including components in blue
# ------------------------------------------------------------------------------
df_sorted <- df |>
  arrange(desc(Percent_of_R2_TMPG)) |>
  mutate(
    Atom = factor(Atom, levels = Atom),
    Highlight = ifelse(str_detect(Atom, "T"), "T-including", "Other")
  )

p_atoms <- ggplot(df_sorted, aes(x = Atom, y = Percent_of_R2_TMPG, fill = Highlight)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = c("T-including" = HILITE_BLUE, "Other" = GREY_OTHER)) +
  labs(
    x = "Non-overlapping variance components",
    y = "Unique contribution to R²\nin full TMPG-model (%)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    axis.title.y = element_text(size = 15, margin = margin(r = 15))  
  )


ggsave(out_atoms,  p_atoms,  width = 8, height = 4, dpi = 300) 
message("Wrote: ", out_atoms)

# ------------------------------------------------------------------------------
# 2) Totals by set (including overlaps), highlight T in blue, others grey
# ------------------------------------------------------------------------------
sets <- c("T","M","P","G")
totals <- lapply(sets, function(s){
  sum(df$Percent_of_R2_TMPG[str_detect(df$Atom, s)], na.rm = TRUE)
})

totals_df <- tibble(
  Set = factor(sets, levels = sets),
  Percent_R2_Total = unlist(totals)
)

# Nicer x-axis labels
label_map <- c(T = "Tetra", M = "Mech", P = "Phylo", G = "Genome")
totals_df$Set_nice <- factor(recode(as.character(totals_df$Set), !!!label_map),
                             levels = label_map)

# Highlight mapping
totals_df <- totals_df |>
  mutate(Highlight = ifelse(Set == "T", "T", "Other"))

p_totals <- ggplot(totals_df, aes(x = Set_nice, y = Percent_R2_Total, fill = Highlight)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = c("T" = HILITE_BLUE, "Other" = GREY_OTHER)) +
  labs(
    x = "Predictors (including overlaps)",
    y = "Total contribution to R²\nin full TMPG-model (%)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_text(size = 15, margin = margin(r = 15)) 
  )


ggsave(out_totals, p_totals, width = 8, height = 4, dpi = 300)
message("Wrote: ", out_totals)

print(totals_df)

