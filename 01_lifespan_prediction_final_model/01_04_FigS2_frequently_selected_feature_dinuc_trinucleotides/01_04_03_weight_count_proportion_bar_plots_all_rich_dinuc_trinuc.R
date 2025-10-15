
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# Tag to keep filenames consistent with the Selection_Frequency filter (> 0.1)
sel_tag <- "sel0.1"

# ---------- Helpers ----------
pretty_trait <- function(x) {
  x |>
    str_replace("^AT_rich$", "AT-rich") |>
    str_replace("^GC_rich$", "GC-rich") |>
    str_replace("^Purine_rich$", "Purine-rich") |>
    str_replace("^Pyrimidine_rich$", "Pyrimidine-rich") |>
    str_replace("^AT_exclusive$", "AT-exclusive") |>
    str_replace("^([ACGT]{2,3})_including$", "\\1-including")
}

# General plotting helper
make_plot_for_traits <- function(props_df, traits_in_order, ncol_facets, outfile_stub) {
  # Keep only traits that actually exist
  traits_present <- traits_in_order[traits_in_order %in% unique(props_df$Trait)]
  if (length(traits_present) == 0) {
    stop("None of the requested traits are present in the props data.")
  }

  # Trait info (for Possible_Tetranucleotide_Number & Expected_Prop)
  trait_info <- props_df |>
    distinct(Trait, Possible_Tetranucleotide_Number, Expected_Prop)

  # Long format for plotting
  props_long <- props_df |>
    filter(Trait %in% traits_present) |>
    mutate(Trait = factor(Trait, levels = traits_present)) |>
    select(Trait, Region, Positive_Prop, Negative_Prop) |>
    pivot_longer(cols = c(Positive_Prop, Negative_Prop),
                 names_to = "Polarity", values_to = "Proportion") |>
    mutate(
      Polarity = ifelse(Polarity == "Positive_Prop", "Positive", "Negative"),
      Region = factor(
        Region,
        levels = c("genome","upstream","exon","intron","downstream"),
        labels = c("Genome","Upstream","Exon","Intron","Downstream")
      )
    ) |>
    left_join(trait_info, by = "Trait")

  # Build label lookup: pretty name + (Possible_Tetranucleotide_Number)
  label_df <- trait_info |>
    filter(Trait %in% traits_present) |>
    mutate(label = paste0(pretty_trait(Trait),
                          " (", Possible_Tetranucleotide_Number, ")")) |>
    select(Trait, label)

  label_lookup <- setNames(label_df$label, label_df$Trait)

  # Build plot
  p <- ggplot(props_long, aes(x = Region, y = Proportion, fill = Polarity)) +
    geom_col(position = "dodge", alpha = 0.9, na.rm = TRUE) +
    facet_wrap(
      ~ Trait,
      ncol = ncol_facets,
      scales = "fixed",
      labeller = labeller(Trait = label_lookup)
    ) +
    scale_fill_manual(values = c("Positive" = "#d73027", "Negative" = "#4575b4")) +
    scale_y_continuous(limits = c(0, 1)) +
    geom_hline(
      aes(yintercept = Expected_Prop),
      data = props_long |> distinct(Trait, Expected_Prop),
      linetype = "dotted", color = "black", size = 1
    ) +
    theme_minimal(base_size = 21) +
    labs(
      x = "Region",
      y = "Proportion of features\n(sel.freq. > 0.1)",
      fill = "Polarity"
    ) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 16),
      strip.text   = element_text(face = "bold", size = 17),
      legend.position = "right",
      axis.title.y = element_text(margin = margin(r = 15))
    )

  # Autosize by grid geometry
  ncols <- ncol_facets
  nrows <- ceiling(length(traits_present) / ncols)
  fig_width  <- max(10, ncols * 3.0)
  fig_height <- max(5,  nrows * 3.0)

  outpath <- paste0("barplot_count_proportions_by_trait_and_region_", sel_tag, "_", outfile_stub, ".png")
  ggsave(outpath, p, width = fig_width, height = fig_height, dpi = 300)
  message("Saved: ", outpath)
  print(p)
}

# ---------- Load props once ----------
props <- read_csv(paste0("count_proportions_by_trait_and_region_", sel_tag, ".csv"),
                  show_col_types = FALSE)

# ---------- Trait sets in desired order ----------
# 1) Richness row (exact order)
richness_traits <- c("AT_exclusive","AT_rich","GC_rich","Purine_rich","Pyrimidine_rich")

# 2) All 16 dinucleotides (8 columns → 2 rows)
dinucs <- c("AA","AC","AG","AT","CA","CC","CG","CT",
            "GA","GC","GG","GT","TA","TC","TG","TT")
dinuc_traits <- paste0(dinucs, "_including")

# 3) All 64 trinucleotides (8 columns → 8 rows)
bases <- c("A","C","G","T")
trinucs <- as.vector(outer(outer(bases, bases, paste0), bases, paste0))
trinuc_traits <- paste0(trinucs, "_including")

# ---------- Make the three PNGs ----------
make_plot_for_traits(props, richness_traits, ncol_facets = 5, outfile_stub = "richness")
make_plot_for_traits(props, dinuc_traits,   ncol_facets = 8, outfile_stub = "dinucleotides_16")
make_plot_for_traits(props, trinuc_traits,  ncol_facets = 8, outfile_stub = "trinucleotides_64")
