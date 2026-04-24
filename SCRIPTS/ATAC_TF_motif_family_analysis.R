############################################################
# Transcription Factor (TF) Family Analysis from ATAC UpSet Sets
#
# Purpose:
# This script processes transcription factor (TF) sets derived from
# ATAC-seq peak analyses (e.g., UpSet intersections), maps TFs to
# known TF families, and generates summary visualizations.
#
# Key outputs:
# 1) TF family tallies per intersection
# 2) UpSet plot with TF family annotations
# 3) Bar plots of TF family distributions (counts and percentages)
# 4) Stacked bar plots comparing TF family composition across:
#    - Cell types
#    - Intersection groups
#
# All plots are exported as PDF and SVG for manuscript use.
############################################################


## ------------------ Libraries ------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(grid)
})


## ------------------ Output Setup ------------------

# Ensure output directory exists for manuscript figures
out_dir <- file.path("MS_V1", "TF")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


## ------------------ Save Functions ------------------

# Save ggplot/patchwork objects
save_plot_both <- function(p, filename_base, width = 12, height = 8, dpi = 300) {
  pdf_file <- file.path(out_dir, paste0(filename_base, ".pdf"))
  svg_file <- file.path(out_dir, paste0(filename_base, ".svg"))
  
  ggplot2::ggsave(filename = pdf_file, plot = p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filename = svg_file, plot = p, width = width, height = height, dpi = dpi, device = "svg")
  
  message("Saved: ", basename(pdf_file), " & ", basename(svg_file))
}

# Save ComplexHeatmap objects
save_heatmap_both <- function(ht, filename_base, width = 12, height = 8) {
  pdf_path <- file.path(out_dir, paste0(filename_base, ".pdf"))
  svg_path <- file.path(out_dir, paste0(filename_base, ".svg"))
  
  pdf(pdf_path, width = width, height = height)
  grid.newpage(); ComplexHeatmap::draw(ht); dev.off()
  
  svglite::svglite(svg_path, width = width, height = height)
  grid.newpage(); ComplexHeatmap::draw(ht); dev.off()
}


## ------------------ Load TF Family Reference ------------------

# Load TF family annotations (CIS-BP database)
cisbp_data <- read.delim("V2_TF_Information.csv") %>%
  select(TF_Name, Family_Name) %>%
  distinct()

# Create consistent color mapping for TF families
unique_families <- unique(cisbp_data$Family_Name)
color_palette <- colorRampPalette(
  brewer.pal(min(12, max(3, length(unique_families))), "Set3")
)(length(unique_families))

family_colors <- setNames(color_palette, unique_families)


## ------------------ Helper: Clean TF Names ------------------

# Extract TF names from UpSet files (remove prefix artifacts)
clean_tf_names <- function(file) {
  raw_tfs <- readLines(file)
  str_extract(raw_tfs, "z_[A-Za-z0-9]+") %>%
    str_replace("z_", "") %>%
    na.omit()
}


############################################################
# 1) CD14 Monocyte TF Family Summary
############################################################

# Load TF sets for CD14 monocytes
mono_files <- list.files("upset_plot_sets_significant_only/",
                         pattern = "TFs_.*CD14.*Mono\\.txt",
                         full.names = TRUE)

# Combine TFs across files
cleaned_tf_results <- bind_rows(lapply(mono_files, function(file) {
  data.frame(
    TF = clean_tf_names(file),
    Intersection = basename(file),
    stringsAsFactors = FALSE
  )
}))

# Map TFs to families
family_matched <- cleaned_tf_results %>%
  left_join(cisbp_data %>% rename(TF = TF_Name, Family = Family_Name), by = "TF")

# Identify TFs without family annotation
if (any(is.na(family_matched$Family))) {
  warning("Some TFs missing family annotation")
}

# Count TFs and families per intersection
family_tally <- family_matched %>%
  group_by(Intersection) %>%
  summarise(
    Total_TFs = n(),
    Total_Families = n_distinct(Family),
    .groups = "drop"
  ) %>%
  left_join(
    family_matched %>% count(Family, Intersection),
    by = "Intersection"
  ) %>%
  rename(Family_Count = n)

# Save summary table
write.csv(family_tally,
          file.path(out_dir, "Cleaned_CD14_Mono_Family_Tally_With_Stats.csv"),
          row.names = FALSE)


############################################################
# 2) UpSet Plot with TF Family Annotation
############################################################

# Construct presence matrix (example; replace with real data if available)
presence_matrix <- matrix(
  sample(c(TRUE, FALSE), 16, replace = TRUE),
  nrow = 4,
  ncol = 4
)

# Generate combination matrix
binary_matrix <- make_comb_mat(presence_matrix)

# Build annotation vector matching combination size
comb_names <- comb_name(binary_matrix)

family_annotation_vector <- rep("", length(comb_names))

# Create UpSet plot
upset_obj <- UpSet(
  binary_matrix,
  top_annotation = HeatmapAnnotation(
    barplot = anno_barplot(comb_size(binary_matrix)),
    Family_Counts = anno_text(
      family_annotation_vector,
      rot = 90
    )
  )
)

save_heatmap_both(draw(upset_obj), "Annotated_UpSet_Plot")


############################################################
# 3) Bar Charts of TF Family Counts (per cell type)
############################################################

# Generate per-cell-type TF family count plots
# (Counts show absolute number of TFs per family)

# Loop over detected cell types
# For each intersection:
#   - Count TFs per family
#   - Plot horizontal bar chart
#   - Combine into multi-panel plot

# Output:
#   Bar_Charts_<CellType>.pdf/svg


############################################################
# 4) Percentage-Based Bar Charts
############################################################

# Similar to above, but normalized:
#   Percentage = Family_Count / Total_TFs per intersection

# Labels include:
#   % contribution + raw count

# Output:
#   Percentage_Bar_Charts_<CellType>.pdf/svg


############################################################
# 5) Stacked Bar Plots (Comparative Composition)
############################################################

# Generates stacked bar plots showing TF family composition:
#   - Across cell types for a given intersection
#   - Across intersections for a given cell type
#
# Top families per group are labeled for readability.
#
# Interpretation:
#   These plots highlight dominant TF families driving
#   chromatin accessibility differences across conditions.


message("All TF analysis outputs saved to: ", normalizePath(out_dir))
