############################################################
# scRNA-seq GSEA Bubble Plot Visualization
#
# Purpose:
# This script visualizes filtered GSEA results generated from the
# pseudobulk delta-log2FC pipeline. It produces bubble plots showing
# pathway enrichment across training conditions.
#
# Plot encoding:
# - X-axis: Training condition (SYKi, MDP, BG)
# - Y-axis: Pathway (Description)
# - Color: Normalized Enrichment Score (NES)
# - Size: Statistical significance (-log10 adjusted p-value)
#
# Filtering:
# - Retains pathways with adjusted p-value < 0.05
# - Removes pathways labeled "TBA"
#
# Output:
# - Separate plots for:
#     1) MDP vs BG
#     2) SYKi only
#     3) All conditions combined
# - Saved as PNG and SVG for manuscript use
############################################################


# ---- Libraries ----
library(ggplot2)
library(dplyr)
library(stringr)


# ---- Paths & Output Dir ----

# Input file containing combined GSEA results
combined_results_path <- file.path("./CURRENT_DATA/", "Filtered_GSEA_Results_v3.csv")

message("\nReading GSEA results from: ", combined_results_path, "\n")

# Output directory for manuscript figures
out_dir <- "./MS_V1/pathway"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ---- Read Data ----
df <- read.csv(combined_results_path, stringsAsFactors = FALSE)


# ---- Cell Type Ordering ----
# Ensures consistent facet order across plots
custom_order <- c(
  "Monocytes", "T cells", "CD4+ T cells", "CD8+ T cells",
  "CD4 TCM", "CD8 TCM", "NK cells", "Dendritic cells", "B cells"
)


# ---- Plotting Function ----
# Generates a faceted bubble plot of pathway enrichment
make_plot <- function(dat) {
  
  # Apply consistent cell type ordering
  dat$CellType <- factor(dat$CellType, levels = custom_order)
  
  ggplot(dat, aes(
    x = Condition,
    y = Description,
    color = NES,
    size = -log10(p.adjust)
  )) +
    geom_point(alpha = 0.8) +
    
    # Diverging color scale centered at NES = 0
    scale_color_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0
    ) +
    
    # Bubble size reflects statistical significance
    scale_size_continuous(range = c(4, 26)) +
    
    labs(
      x = "Condition",
      y = "Pathway",
      color = "NES",
      size = "-log10(p.adjust)",
      title = "Pathway Enrichment Analysis"
    ) +
    
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16),
      axis.text.y = element_text(size = 18),
      strip.text  = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      panel.spacing = unit(3, "lines")
    ) +
    
    # Wrap long pathway names for readability
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    
    # Facet by cell type
    facet_wrap(~ CellType, scales = "free_y", ncol = 3) +
    
    coord_cartesian(clip = "off")
}


# ---- Save Function ----
# Saves plots in both PNG and SVG formats
save_both <- function(plot, basename,
                      width = 24, height = 34, dpi = 300,
                      units = "in", out_dir = getwd()) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  svg_path <- file.path(out_dir, paste0(basename, ".svg"))
  
  ggplot2::ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi
  )
  
  ggplot2::ggsave(
    filename = svg_path,
    plot = plot,
    width = width,
    height = height,
    units = units,
    device = "svg"
  )
  
  message("Saved: ", basename, " (.png & .svg) -> ", out_dir)
}


############################################################
# 1) MDP + BG Comparison
############################################################

# Identify significant pathways in MDP and BG
sig_ids_mdp_bg <- df %>%
  filter(Condition %in% c("MDP", "BG")) %>%
  filter(p.adjust < 0.05) %>%
  filter(!str_detect(ID, "TBA")) %>%
  distinct(ID) %>%
  pull(ID)

# Subset dataset to these pathways
df_mdp_bg <- df %>%
  filter(ID %in% sig_ids_mdp_bg, Condition %in% c("MDP", "BG"))

# Generate and save plot
plot_mdp_bg <- make_plot(df_mdp_bg)
save_both(plot_mdp_bg, "v03_Pathway_Bubble_MDP_BG", out_dir = out_dir)


############################################################
# 2) SYKi Only
############################################################

sig_ids_syki <- df %>%
  filter(Condition == "SYKi") %>%
  filter(p.adjust < 0.05) %>%
  filter(!str_detect(ID, "TBA")) %>%
  distinct(ID) %>%
  pull(ID)

df_syki <- df %>%
  filter(ID %in% sig_ids_syki, Condition == "SYKi")

plot_syki <- make_plot(df_syki)
save_both(plot_syki, "v03_Pathway_Bubble_SYKi", out_dir = out_dir)


############################################################
# 3) All Conditions Combined
############################################################

sig_ids_all <- df %>%
  filter(p.adjust < 0.05) %>%
  filter(!str_detect(ID, "TBA")) %>%
  distinct(ID) %>%
  pull(ID)

df_all <- df %>%
  filter(ID %in% sig_ids_all)

plot_all <- make_plot(df_all)
save_both(plot_all, "v03_Pathway_Bubble_ALL", out_dir = out_dir)


message("All plots saved (PNG + SVG) to: ", out_dir, "\n")
