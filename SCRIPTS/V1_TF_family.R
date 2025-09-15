###############################################
# ONE CONTINUOUS UPDATED COMPLETE R SCRIPT
# - Saves all plots to MS_V1/TF/ as both PDF & SVG
# - Robust fix for ComplexHeatmap anno_text length mismatch
###############################################

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

## ------------------ Output helpers ------------------
# Ensure svglite is available for SVG saving
if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}
out_dir <- file.path("MS_V1", "TF")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save ggplot/patchwork objects to PDF + SVG
save_plot_both <- function(p, filename_base, width = 12, height = 8, dpi = 300) {
  pdf_file <- file.path(out_dir, paste0(filename_base, ".pdf"))
  svg_file <- file.path(out_dir, paste0(filename_base, ".svg"))
  ggplot2::ggsave(filename = pdf_file, plot = p, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(filename = svg_file, plot = p, width = width, height = height, dpi = dpi, device = "svg")
  message("Saved: ", basename(pdf_file), " & ", basename(svg_file))
}

# Save ComplexHeatmap / grid object returned by draw() to PDF + SVG
save_heatmap_both <- function(ht, filename_base, width = 12, height = 8) {
  pdf_path <- file.path(out_dir, paste0(filename_base, ".pdf"))
  svg_path <- file.path(out_dir, paste0(filename_base, ".svg"))
  
  grDevices::pdf(pdf_path, width = width, height = height)
  grid::grid.newpage(); ComplexHeatmap::draw(ht); grDevices::dev.off()
  
  svglite::svglite(svg_path, width = width, height = height)
  grid::grid.newpage(); ComplexHeatmap::draw(ht); grDevices::dev.off()
  
  message("Saved: ", basename(pdf_path), " & ", basename(svg_path))
}

# CSV path helper (optional)
csv_out <- function(name) file.path(out_dir, name)

## ------------------ Common inputs ------------------
upset_dir <- "upset_plot_sets_significant_only/"

# Load cisbp_data (tries CSV first, then TXT)
cisbp_path_csv <- "V2_TF_Information.csv"
cisbp_path_txt <- "TF_Information.txt"

if (file.exists(cisbp_path_csv)) {
  cisbp_data <- read.delim(cisbp_path_csv, header = TRUE, sep = ",")
} else if (file.exists(cisbp_path_txt)) {
  cisbp_data <- read.delim(cisbp_path_txt, header = TRUE, sep = "\t")
} else {
  stop("Could not find V2_TF_Information.csv or TF_Information.txt for TF family mapping.")
}
cisbp_data <- cisbp_data %>% dplyr::select(TF_Name, Family_Name) %>% distinct()

# Global TF-family color mapping
unique_families <- unique(cisbp_data$Family_Name)
color_palette <- colorRampPalette(brewer.pal(min(12, max(3, length(unique_families))), "Set3"))(length(unique_families))
family_colors <- setNames(color_palette, unique_families)

# Utility: Read TF names from a file and clean
clean_tf_names <- function(file) {
  raw_tfs <- readLines(file)
  clean_tfs <- str_extract(raw_tfs, "z_[A-Za-z0-9]+") %>%
    str_replace("z_", "") %>%
    na.omit()
  return(clean_tfs)
}

## ===================================================
## 1) CD14 Mono: Parse, tally, validate, save CSV
## ===================================================

# List significant files for CD14 Mono
mono_files <- list.files(upset_dir, pattern = "TFs_.*CD14.*Mono\\.txt", full.names = TRUE)

# Parse & clean
cleaned_tf_data <- lapply(mono_files, function(file) {
  data.frame(
    TF = clean_tf_names(file),
    Intersection = basename(file),
    stringsAsFactors = FALSE
  )
})
cleaned_tf_results <- bind_rows(cleaned_tf_data)

# Join with cisbp families
family_matched <- cleaned_tf_results %>%
  left_join(
    cisbp_data %>% dplyr::rename(TF = TF_Name, Family = Family_Name),
    by = "TF"
  )

# Warn on unmatched
unmatched_tfs <- family_matched %>% filter(is.na(Family))
if (nrow(unmatched_tfs) > 0) {
  warning("Some TFs do not have families assigned. Review 'unmatched_tfs' object.")
  print(unmatched_tfs)
}

# Family tally per intersection
family_tally <- family_matched %>%
  group_by(Intersection) %>%
  summarise(
    Total_TFs = n(),
    Total_Families = n_distinct(Family),
    .groups = "drop"
  ) %>%
  left_join(
    family_matched %>%
      count(Family, Intersection, sort = TRUE),
    by = "Intersection"
  ) %>%
  dplyr::rename(Family_Count = n)

# Validate TF counts
total_tf_count <- nrow(cleaned_tf_results)
total_family_tally <- sum(family_tally$Family_Count, na.rm = TRUE)
if (total_tf_count != total_family_tally) {
  stop(paste("Mismatch in TF count and family tally:",
             "Total TFs =", total_tf_count,
             "Total Family Tally =", total_family_tally))
} else {
  message("TF count matches family tally. Validation passed.")
}

# Save CSV to MS_V1/TF/
write.csv(family_tally, csv_out("Cleaned_CD14_Mono_Family_Tally_With_Stats.csv"), row.names = FALSE)

## ===================================================
## 2) ComplexHeatmap UpSet example with family text
##    (FIX: make label vector exactly length of combinations)
## ===================================================

# Example presence matrix (replace with real presence matrix if available)
presence_matrix <- matrix(
  data = sample(c(TRUE, FALSE), 16, replace = TRUE),
  nrow = 4,
  ncol = 4,
  dimnames = list(
    c("TFs_000_CD14 Mono.txt", "TFs_001_CD14 Mono.txt", "TFs_010_CD14 Mono.txt", "TFs_011_CD14 Mono.txt"),
    c("C2H2 ZF", "bZIP", "IRF", "bHLH")
  )
)

# Build combination matrix + names
binary_matrix <- make_comb_mat(presence_matrix)
comb_names <- comb_name(binary_matrix)                # character vector of combination names
n_combs <- length(comb_names)

# Build a named lookup for annotations from family_tally
# (Intersections in presence_matrix rownames should match filenames in family_tally$Intersection)
ann_lookup <- family_tally %>%
  group_by(Intersection) %>%
  summarise(Family_Annotation = paste(Family, "(", Family_Count, ")", collapse = "; "), .groups = "drop") %>%
  { setNames(.$Family_Annotation, .$Intersection) }

# Create a label for every combination column, in order; fill missing with ""
# NOTE: In UpSet, combo columns correspond to set intersections, while our
#       annotations are per-file/row intersections; if they don't align, we just show blanks.
family_annotation_vector <- rep("", n_combs)
names(family_annotation_vector) <- comb_names
# If your comb_names are actual filenames, this will fill; otherwise remains ""
has_match <- intersect(names(ann_lookup), comb_names)
family_annotation_vector[has_match] <- ann_lookup[has_match]

# Construct UpSet with length-matching anno_text
upset_obj <- UpSet(
  binary_matrix,
  set_order = colnames(presence_matrix),
  top_annotation = HeatmapAnnotation(
    barplot = anno_barplot(comb_size(binary_matrix), gp = gpar(fill = "steelblue")),
    Family_Counts = anno_text(
      family_annotation_vector,     # exactly n_combs long
      rot = 90,
      just = "center",
      gp = gpar(fontsize = 8, col = "darkred")
    )
  ),
  row_names_gp = gpar(fontsize = 10),
  column_title = "UpSet Plot with Family Annotations",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

upset_drawn <- draw(upset_obj)
save_heatmap_both(upset_drawn, "Annotated_UpSet_Plot_with_Family_Counts", width = 12, height = 8)

## ===================================================
## 3) Per-cell-type bar charts (counts)
## ===================================================

all_files <- list.files(upset_dir, pattern = "TFs_.*\\.txt", full.names = TRUE)
cell_types <- unique(str_extract(all_files, "(?<=TFs_\\d{3}_).*(?=\\.txt)"))

for (cell_type in cell_types) {
  celltype_files <- all_files[str_detect(all_files, fixed(cell_type))]
  cleaned_tf_data_ct <- lapply(celltype_files, function(file) {
    data.frame(
      TF = clean_tf_names(file),
      Intersection = basename(file),
      stringsAsFactors = FALSE
    )
  })
  cleaned_tf_results_ct <- bind_rows(cleaned_tf_data_ct)
  
  if (nrow(cleaned_tf_results_ct) == 0) next
  
  family_matched_ct <- cleaned_tf_results_ct %>%
    left_join(
      cisbp_data %>% dplyr::rename(TF = TF_Name, Family = Family_Name),
      by = "TF"
    )
  
  family_tally_ct <- family_matched_ct %>%
    dplyr::group_by(Intersection) %>%
    dplyr::summarise(
      Total_TFs = n(),
      Total_Families = n_distinct(Family),
      .groups = "drop"
    ) %>%
    left_join(
      family_matched_ct %>% dplyr::count(Family, Intersection, sort = TRUE),
      by = "Intersection"
    ) %>%
    dplyr::rename(Family_Count = n)
  
  bar_charts <- lapply(unique(family_tally_ct$Intersection), function(intersection) {
    data <- family_tally_ct %>% filter(Intersection == intersection)
    total_count <- sum(data$Family_Count, na.rm = TRUE)
    ggplot(data, aes(x = reorder(Family, -Family_Count), y = Family_Count, fill = Family)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Family_Count), vjust = -0.5, size = 3) +
      coord_flip() +
      labs(
        title = paste(intersection, "- Total:", total_count),
        x = "Family", y = "Count"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  combined_plot <- wrap_plots(bar_charts, ncol = 2)
  save_plot_both(combined_plot, paste0("Bar_Charts_", cell_type), width = 16, height = 12)
  message("Saved combined bar charts for cell type: ", cell_type)
}

## ===================================================
## 4) Per-cell-type bar charts (percentages + counts in labels)
## ===================================================

for (cell_type in cell_types) {
  celltype_files <- all_files[str_detect(all_files, fixed(cell_type))]
  cleaned_tf_data_ct <- lapply(celltype_files, function(file) {
    data.frame(
      TF = clean_tf_names(file),
      Intersection = basename(file),
      stringsAsFactors = FALSE
    )
  })
  cleaned_tf_results_ct <- bind_rows(cleaned_tf_data_ct)
  if (nrow(cleaned_tf_results_ct) == 0) next
  
  family_matched_ct <- cleaned_tf_results_ct %>%
    left_join(
      cisbp_data %>% dplyr::rename(TF = TF_Name, Family = Family_Name),
      by = "TF"
    )
  
  family_tally_ct <- family_matched_ct %>%
    group_by(Intersection) %>%
    dplyr::summarise(
      Total_TFs = n(),
      Total_Families = n_distinct(Family),
      .groups = "drop"
    ) %>%
    left_join(
      family_matched_ct %>% dplyr::count(Family, Intersection, sort = TRUE),
      by = "Intersection"
    ) %>% dplyr::rename(Family_Count = n) %>%
    group_by(Intersection) %>%
    mutate(Total_TFs_Per_Intersection = sum(Family_Count, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      Percentage = 100 * Family_Count / Total_TFs_Per_Intersection,
      Label = paste0(sprintf("%.1f", Percentage), "% (", Family_Count, ")")
    )
  
  bar_charts <- lapply(unique(family_tally_ct$Intersection), function(intersection) {
    data <- family_tally_ct %>% filter(Intersection == intersection)
    total_count <- unique(data$Total_TFs_Per_Intersection)
    ggplot(data, aes(x = reorder(Family, -Percentage), y = Percentage, fill = Family)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = Label), vjust = -0.5, size = 3) +
      coord_flip() +
      labs(
        title = paste(intersection, "- Total TFs:", total_count),
        x = "Family", y = "Percentage of Total TFs (%)"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  })
  
  combined_plot <- wrap_plots(bar_charts, ncol = 2)
  save_plot_both(combined_plot, paste0("Percentage_Bar_Charts_with_Counts_", cell_type), width = 16, height = 16)
  message("Saved percentage-based combined bar charts for cell type: ", cell_type)
}

## ===================================================
## 5) Stacked bars by INTERSECTION PREFIX (v3_Stacked_*)
## ===================================================

# Re-list files in case earlier blocks changed environment
all_files <- list.files(upset_dir, pattern = "TFs_.*\\.txt", full.names = TRUE)
intersection_prefixes <- unique(str_extract(all_files, "TFs_\\d{3}"))

# (a) v3_Stacked_Bar_Chart_<prefix>.*
for (prefix in intersection_prefixes) {
  matching_files <- all_files[str_detect(all_files, prefix)]
  cleaned_tf_data_px <- lapply(matching_files, function(file) {
    data.frame(
      TF = clean_tf_names(file),
      Intersection = basename(file),
      CellType = str_extract(basename(file), "(?<=TFs_\\d{3}_).*(?=\\.txt)"),
      stringsAsFactors = FALSE
    )
  })
  combined_data <- bind_rows(cleaned_tf_data_px)
  if (nrow(combined_data) == 0) { message("No data for prefix: ", prefix); next }
  
  family_matched <- combined_data %>%
    left_join(cisbp_data %>% dplyr::rename(TF = TF_Name, Family = Family_Name), by = "TF")
  
  family_tally <- family_matched %>%
    group_by(CellType) %>%
    dplyr::summarise(Total_TFs = n(), .groups = "drop") %>%
    left_join(
      family_matched %>% dplyr::count(Family, CellType, sort = TRUE),
      by = "CellType"
    ) %>% dplyr::rename(Family_Count = n) %>%
    group_by(CellType) %>%
    mutate(
      Total_TFs_Per_CellType = sum(Family_Count, na.rm = TRUE),
      Percentage = 100 * Family_Count / Total_TFs_Per_CellType
    ) %>% ungroup() %>%
    group_by(CellType) %>%
    mutate(
      RankWithinBar = rank(-Percentage, ties.method = "first"),
      Label = ifelse(RankWithinBar <= 3 & !is.na(Family),
                     paste0(Family, "\n", Family_Count, " (", sprintf("%.1f", Percentage), "%)"),
                     NA)
    ) %>% ungroup() %>%
    mutate(CellTypeLabel = paste0(CellType, " (", Total_TFs_Per_CellType, ")"))
  
  p_v3 <- ggplot(family_tally, aes(x = reorder(CellTypeLabel, -Total_TFs_Per_CellType), y = Percentage, fill = Family)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Label, group = Family),
              position = position_stack(vjust = 0.5), size = 3, na.rm = TRUE) +
    scale_fill_manual(values = family_colors) +
    labs(title = paste("Intersection:", prefix),
         x = "Cell Type (Total TFs)", y = "Percentage of Total TFs (%)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  save_plot_both(p_v3, paste0("v3_Stacked_Bar_Chart_", prefix), width = 15, height = 8)
}

# (b) v3C_Stacked_Bar_Chart_<cell_type>.*
for (cell_type in cell_types) {
  celltype_files <- all_files[str_detect(all_files, fixed(cell_type))]
  cleaned_tf_data_ct <- lapply(celltype_files, function(file) {
    data.frame(
      TF = clean_tf_names(file),
      Intersection = str_extract(basename(file), "TFs_\\d{3}"),
      CellType = cell_type,
      stringsAsFactors = FALSE
    )
  })
  combined_data <- bind_rows(cleaned_tf_data_ct)
  if (nrow(combined_data) == 0) { next }
  
  family_matched <- combined_data %>%
    left_join(cisbp_data %>% dplyr::rename(TF = TF_Name, Family = Family_Name), by = "TF")
  
  family_tally <- family_matched %>%
    dplyr::group_by(Intersection, Family) %>%
    dplyr::summarise(Family_Count = n(), .groups = "drop") %>%
    dplyr::group_by(Intersection) %>%
    dplyr::mutate(
      Total_TFs_Per_Intersection = sum(Family_Count, na.rm = TRUE),
      Percentage = 100 * Family_Count / Total_TFs_Per_Intersection,
      Valid = !is.na(Family),
      AdjustedRank = rank(-Percentage * Valid, ties.method = "first"),
      Label = ifelse(AdjustedRank <= 3 & Valid,
                     paste0(Family, "\n", Family_Count, " (", sprintf("%.1f", Percentage), "%)"),
                     NA)
    ) %>% ungroup() %>%
    mutate(IntersectionLabel = paste0(Intersection, " (", Total_TFs_Per_Intersection, ")"))
  
  p_v3c <- ggplot(family_tally, aes(x = reorder(IntersectionLabel, -Total_TFs_Per_Intersection), y = Percentage, fill = Family)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Label, group = Family),
              position = position_stack(vjust = 0.5), size = 3, na.rm = TRUE) +
    scale_fill_manual(values = family_colors) +
    labs(title = paste("Cell Type:", cell_type),
         x = "Intersection (Total TFs)", y = "Percentage of Total TFs (%)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  save_plot_both(p_v3c, paste0("v3C_Stacked_Bar_Chart_", cell_type), width = 15, height = 8)
}

# (c) v3B_Stacked_Bar_Chart_Updated_<prefix>.* and v3C_Stacked_Bar_Chart_Updated_<prefix>.*
for (prefix in intersection_prefixes) {
  matching_files <- all_files[str_detect(all_files, prefix)]
  cleaned_tf_data_px <- lapply(matching_files, function(file) {
    data.frame(
      TF = clean_tf_names(file),
      Intersection = basename(file),
      CellType = str_extract(basename(file), "(?<=TFs_\\d{3}_).*(?=\\.txt)"),
      stringsAsFactors = FALSE
    )
  })
  combined_data <- bind_rows(cleaned_tf_data_px)
  if (nrow(combined_data) == 0) { next }
  
  family_matched <- combined_data %>%
    left_join(cisbp_data %>% dplyr::rename(TF = TF_Name, Family = Family_Name), by = "TF")
  
  family_tally <- family_matched %>%
    dplyr::group_by(CellType, Intersection, Family) %>%
    dplyr::summarise(Family_Count = n(), .groups = "drop") %>%
    dplyr::group_by(CellType, Intersection) %>%
    dplyr::mutate(
      Total_TFs_Per_CellType = sum(Family_Count, na.rm = TRUE),
      Percentage = 100 * Family_Count / Total_TFs_Per_CellType
    ) %>% ungroup() %>%
    dplyr::group_by(CellType, Intersection) %>%
    dplyr::mutate(
      Valid = !is.na(Family),
      AdjustedRank = rank(-Percentage * Valid, ties.method = "first"),
      Label = ifelse(AdjustedRank <= 3 & Valid,
                     paste0(Family, "\n", Family_Count, " (", sprintf("%.1f", Percentage), "%)"),
                     NA)
    ) %>% ungroup() %>%
    dplyr::mutate(CellTypeLabel = paste0(CellType, " (", Total_TFs_Per_CellType, ")"))
  
  plot2 <- ggplot(family_tally, aes(x = reorder(CellTypeLabel, -Total_TFs_Per_CellType), y = Percentage, fill = Family)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Label, group = Family),
              position = position_stack(vjust = 0.5), size = 3, na.rm = TRUE) +
    scale_fill_manual(values = family_colors) +
    labs(title = paste("Intersection:", prefix),
         x = "Cell Type (Total TFs)", y = "Percentage of Total TFs (%)") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  save_plot_both(plot2, paste0("v3B_Stacked_Bar_Chart_Updated_", prefix), width = 12, height = 8)
  save_plot_both(plot2, paste0("v3C_Stacked_Bar_Chart_Updated_", prefix), width = 12, height = 8)
}

message("All outputs saved to: ", normalizePath(out_dir))
