############################################################
# Pseudobulk scRNA-seq Delta-of-Delta GSEA Pipeline
#
# Purpose:
# This script compares the PMA-induced transcriptional response
# across training conditions: Untrained, SYKi, MDP, and BG.
#
# Main idea:
# 1. For each condition, perform Treated vs Untreated DE analysis.
# 2. Use Untrained as the baseline PMA response.
# 3. For each trained condition, calculate:
#
#    delta_log2FC =
#      log2FC(Treated vs Untreated in trained condition)
#      - log2FC(Treated vs Untreated in Untrained)
#
# 4. Rank genes by delta_log2FC.
# 5. Run GSEA to identify pathways enhanced or reduced by training.
#
# Interpretation:
# Positive delta_log2FC means the trained condition has a stronger
# PMA response than Untrained.
#
# Negative delta_log2FC means the trained condition has a weaker
# PMA response than Untrained.
############################################################


############################################################
# Load libraries
############################################################

library(Seurat)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(tibble)
library(tidyverse)
library(MAST, lib = "./lib")
library(progress, lib = "./lib")


############################################################
# Load input data
############################################################

# Load Seurat object.
obj <- readRDS("./418_OUTPUT_ALL_5000CUT/data2_ref.rds")

# Create timestamped output directory for reproducibility.
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(timestamp)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load custom pathway gene sets.
gmt_file <- "BTM_for_GSEA_20131008.gmt"
custom_gmt <- read.gmt(gmt_file)

# Convert RNA assay to standard Assay object for compatibility.
obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay")


############################################################
# Define experimental design
############################################################

# Odd-numbered lanes are untreated; even-numbered lanes are PMA-treated.
untreated_lanes <- c(1, 3, 5, 7, 33, 35, 37, 39)
treated_lanes   <- c(2, 4, 6, 8, 34, 36, 38, 40)

# Lane groups define training conditions.
lane_groups <- list(
  "untrained" = c(1, 2, 33, 34),
  "SYKi"      = c(3, 4, 35, 36),
  "MDP"       = c(5, 6, 37, 38),
  "BG"        = c(7, 8, 39, 40)
)

# Add treatment label to metadata.
obj$Treatment <- ifelse(obj$Lane %in% untreated_lanes, "Untreated", "Treated")


############################################################
# Generate Untrained baseline DE results
############################################################

# The Untrained condition is used as the baseline PMA response.
untrained_lanes <- lane_groups[["untrained"]]

obj_untrained <- subset(
  obj,
  cells = WhichCells(obj, expression = Lane %in% untrained_lanes)
)

# Aggregate cells into pseudobulk profiles.
# Grouping by sCALL and Lane preserves subject/lane-level structure.
pseudo_bulk_untrained <- AggregateExpression(
  obj_untrained,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("Treatment", "monaco_ann1", "sCALL", "Lane")
)

cell_types <- unique(pseudo_bulk_untrained$monaco_ann1)
untrained_degs <- list()

for (cell_type in cell_types) {
  
  if (is.na(cell_type)) next
  
  subset_cell_type <- subset(
    pseudo_bulk_untrained,
    subset = monaco_ann1 == cell_type
  )
  
  treated_cells <- WhichCells(
    subset_cell_type,
    expression = Treatment == "Treated"
  )
  
  untreated_cells <- WhichCells(
    subset_cell_type,
    expression = Treatment == "Untreated"
  )
  
  # Skip cell types without enough pseudobulk samples in both groups.
  if (length(treated_cells) <= 1 || length(untreated_cells) <= 1) {
    message(
      paste(
        "Skipping cell type", cell_type,
        "in Untrained due to insufficient Treated or Untreated samples"
      )
    )
    next
  }
  
  Idents(subset_cell_type) <- "Treatment"
  
  de_result <- FindMarkers(
    subset_cell_type,
    ident.1 = "Treated",
    ident.2 = "Untreated",
    min.pct = 0.25,
    test.use = "DESeq2"
  )
  
  de_result <- de_result %>%
    rownames_to_column(var = "gene")
  
  untrained_degs[[cell_type]] <- de_result
  
  file_name <- paste0(gsub(" ", "_", cell_type), "_deg_untrained.csv")
  file_name <- gsub("/", "_", file_name)
  file_path <- file.path(output_dir, file_name)
  
  write.csv(de_result, file_path, row.names = FALSE)
  
  cat("Stored Untrained DEGs for cell type:", cell_type, "\n")
}


############################################################
# Compare trained conditions against Untrained baseline
############################################################

combined_gsea_results <- list()

for (condition in setdiff(names(lane_groups), "untrained")) {
  
  lanes <- lane_groups[[condition]]
  
  subset_obj <- subset(
    obj,
    cells = WhichCells(obj, expression = Lane %in% lanes)
  )
  
  pseudo_bulk_trained <- AggregateExpression(
    subset_obj,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = c("Treatment", "monaco_ann1", "sCALL", "Lane")
  )
  
  cell_types <- unique(pseudo_bulk_trained$monaco_ann1)
  
  for (cell_type in cell_types) {
    
    if (is.na(cell_type)) next
    
    # Only compare cell types that were also detected in Untrained.
    if (!(cell_type %in% names(untrained_degs))) next
    
    subset_cell_type <- subset(
      pseudo_bulk_trained,
      subset = monaco_ann1 == cell_type
    )
    
    treated_cells <- WhichCells(
      subset_cell_type,
      expression = Treatment == "Treated"
    )
    
    untreated_cells <- WhichCells(
      subset_cell_type,
      expression = Treatment == "Untreated"
    )
    
    if (length(treated_cells) <= 1 || length(untreated_cells) <= 1) {
      message(
        paste(
          "Skipping cell type", cell_type,
          "in condition", condition,
          "due to insufficient Treated or Untreated samples"
        )
      )
      next
    }
    
    Idents(subset_cell_type) <- "Treatment"
    
    de_result_trained <- FindMarkers(
      subset_cell_type,
      ident.1 = "Treated",
      ident.2 = "Untreated",
      min.pct = 0.25,
      test.use = "DESeq2"
    )
    
    de_result_trained <- de_result_trained %>%
      rownames_to_column(var = "gene")
    
    file_name <- paste0(
      gsub(" ", "_", cell_type),
      "_deg_",
      condition,
      ".csv"
    )
    
    file_name <- gsub("/", "_", file_name)
    file_path <- file.path(output_dir, file_name)
    
    write.csv(de_result_trained, file_path, row.names = FALSE)
    
    
    ########################################################
    # Delta-of-delta calculation
    ########################################################
    
    de_result_untrained <- untrained_degs[[cell_type]]
    
    overlapping_genes <- intersect(
      de_result_untrained$gene,
      de_result_trained$gene
    )
    
    if (length(overlapping_genes) == 0) next
    
    delta_result <- de_result_trained %>%
      dplyr::filter(gene %in% overlapping_genes) %>%
      inner_join(
        de_result_untrained %>%
          dplyr::filter(gene %in% overlapping_genes),
        by = "gene",
        suffix = c(".trained", ".untrained")
      ) %>%
      mutate(
        delta_log2FC = avg_log2FC.trained - avg_log2FC.untrained
      ) %>%
      arrange(desc(delta_log2FC)) %>%
      dplyr::select(gene, delta_log2FC)
    
    ranked_genes <- setNames(
      delta_result$delta_log2FC,
      delta_result$gene
    )
    
    
    ########################################################
    # GSEA using delta-ranked genes
    ########################################################
    
    gsea_result <- tryCatch(
      {
        GSEA(
          geneList = ranked_genes,
          TERM2GENE = custom_gmt,
          pAdjustMethod = "BH",
          pvalueCutoff = 2,
          exponent = 1,
          minGSSize = 10,
          maxGSSize = 500
        )
      },
      error = function(e) {
        message(
          paste(
            "Error in GSEA for condition:", condition,
            "and cell type:", cell_type,
            "-", e$message
          )
        )
        NULL
      }
    )
    
    if (!is.null(gsea_result) && length(gsea_result@result$ID) > 0) {
      
      gsea_df <- as.data.frame(gsea_result)
      
      file_name_gsea <- paste0(
        condition,
        "_",
        gsub(" ", "_", cell_type),
        "_gsea.csv"
      )
      
      file_name_gsea <- gsub("/", "_", file_name_gsea)
      file_path_gsea <- file.path(output_dir, file_name_gsea)
      
      write.csv(gsea_df, file_path_gsea, row.names = FALSE)
      
      gsea_df$Condition <- condition
      gsea_df$CellType <- cell_type
      
      combined_gsea_results[[paste(condition, cell_type, sep = "_")]] <- gsea_df
      
      cat(
        "Stored GSEA for condition:",
        condition,
        "and cell type:",
        cell_type,
        "\n"
      )
      
    } else {
      message(
        paste(
          "No enriched terms for condition:",
          condition,
          "and cell type:",
          cell_type
        )
      )
    }
  }
}


############################################################
# Save combined GSEA results
############################################################

all_gsea_results <- bind_rows(combined_gsea_results)

combined_results_path <- file.path(
  output_dir,
  "GSEA_Pathway_enrichment_unfiltered.csv"
)

write.csv(all_gsea_results, combined_results_path, row.names = FALSE)

cat("\nCombined GSEA results saved to", combined_results_path, "\n")
