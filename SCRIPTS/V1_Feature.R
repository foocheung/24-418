#!/usr/bin/env Rscript
# Plot Seurat FeaturePlot-style condition panels + violin summaries.
# Now supports BOTH:
#   1) Module features (ModuleScore* columns in meta.data)
#   2) Gene features (from RNA assay, e.g., GZMB, IFNG), saved as PDF + SVG
#
# Usage examples:
#   Rscript V1_Feature.R  obj_mod.rds
#   Rscript V1_Feature.R  obj_mod.rds "ModuleScore_Inflam,ModuleScore_IFN"
#   Rscript V1_Feature.R  obj_mod.rds modules.txt
#   Rscript V1_Feature.R  obj_mod.rds "ModuleScore_IFN" modules.txt "GZMB,IFNG"
#   Rscript V1_Feature.R  obj_mod.rds "" "" genes.txt
#
# Notes:
# - If the 2nd arg looks like a file, it's treated as the modules file.
# - Genes can be passed via 4th (comma-separated) and/or 5th (file) args.
# - Saves to ./MS_V1/feature_plot as PDF + SVG. 
#Rscript V1_Feature.R  obj_mod.rds module.txt
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(fs)
  library(rlang)
  library(tibble)
})

# ------------------ CLI ------------------
args <- commandArgs(trailingOnly = TRUE)
rds_path      <- if (length(args) >= 1) args[1] else "obj_mod.rds"
modules_arg   <- if (length(args) >= 2) args[2] else ""  # comma-separated module names or regex
modules_file  <- if (length(args) >= 3) args[3] else ""  # file listing modules/regex (one per line)
genes_arg     <- if (length(args) >= 4) args[4] else ""  # comma-separated gene symbols
genes_file    <- if (length(args) >= 5) args[5] else ""  # file listing gene symbols (one per line)

# If 2nd arg is actually a file path, treat as modules_file
if (nzchar(modules_arg) && fs::file_exists(modules_arg) && !fs::dir_exists(modules_arg)) {
  modules_file <- modules_arg
  modules_arg  <- ""
}

# ------------------ CONFIG ------------------
OUTDIR <- "MS_V1/feature_plot"
fs::dir_create(OUTDIR)

# Label controls for panel labels
ADD_LABELS <- TRUE
MAX_LABELS <- 120L

# Violin median text
SHOW_MEDIAN_TEXT <- TRUE
MEDIAN_DIGITS    <- 2

# Violin fill palette
FILL_VALS <- c(
  "untrained_Untreated" = "lightblue", "SYKi_Untreated" = "blue",
  "MDP_Untreated"       = "purple",    "BG_Untreated"   = "pink",
  "untrained_Treated"   = "lightblue", "SYKi_Treated"   = "blue",
  "MDP_Treated"         = "purple",    "BG_Treated"     = "pink"
)

# Gene defaults if none are provided (only plotted if present in RNA)
DEFAULT_GENES <- c("GZMB", "IFNG")

# Whether to annotate gene violins using NON-ZERO values for dashed reference medians
GENE_VIOLIN_NONZERO_ONLY <- TRUE

# ------------------ HELPERS ------------------
parse_cli_list <- function(s) {
  if (!nzchar(s)) character(0) else trimws(unlist(strsplit(s, ",")))
}
parse_file_list <- function(path) {
  if (!nzchar(path)) return(character(0))
  if (!fs::file_exists(path)) stop("List file not found: ", path)
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  x <- x[!grepl("^\\s*#", x)]
  unique(x)
}
pick_with_regex <- function(all_names, queries) {
  if (length(queries) == 0) return(all_names)
  exact <- intersect(all_names, queries)
  left  <- setdiff(queries, exact)
  regex_hits <- if (length(left)) grep(paste(left, collapse="|"), all_names, ignore.case = TRUE, value = TRUE) else character(0)
  unique(c(exact, regex_hits))
}
# save_pdf_svg <- function(plot_obj, base_name, width = 12, height = 12) {
#   pdf_path <- fs::path(OUTDIR, paste0(base_name, ".pdf"))
#   svg_path <- fs::path(OUTDIR, paste0(base_name, ".svg"))
#   ggsave(filename = pdf_path, plot = plot_obj, width = width, height = height, limitsize = FALSE)
#   ggsave(filename = svg_path, plot = plot_obj, width = width, height = height, device = "svg", limitsize = FALSE)
#   message("Saved: ", basename(pdf_path), " and ", basename(svg_path))
# }


save_pdf_svg <- function(plot_obj,
                         base_name,
                         width = 12,
                         height = 12,
                         units = "in",
                         out_dir = NULL,
                         limitsize = FALSE,
                         use_svglite = TRUE) {
  # Resolve output directory
  if (is.null(out_dir)) {
    out_dir <- if (exists("OUTDIR", inherits = TRUE)) get("OUTDIR", inherits = TRUE) else getwd()
  }
  if (requireNamespace("fs", quietly = TRUE)) {
    fs::dir_create(out_dir)
    pdf_path <- fs::path(out_dir, paste0(base_name, ".pdf"))
    svg_path <- fs::path(out_dir, paste0(base_name, ".svg"))
  } else {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_path <- file.path(out_dir, paste0(base_name, ".pdf"))
    svg_path <- file.path(out_dir, paste0(base_name, ".svg"))
  }
  
  # Always call ggplot2's ggsave and pass a single units string
  ggplot2::ggsave(filename = pdf_path,
                  plot = plot_obj,
                  width = width, height = height, units = units,
                  limitsize = limitsize)
  
  # SVG device: prefer svglite if available
  if (use_svglite && requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(filename = svg_path,
                    plot = plot_obj,
                    width = width, height = height, units = units,
                    device = svglite::svglite,
                    limitsize = limitsize)
  } else {
    ggplot2::ggsave(filename = svg_path,
                    plot = plot_obj,
                    width = width, height = height, units = units,
                    device = "svg",
                    limitsize = limitsize)
  }
  
  message("Saved: ", basename(pdf_path), " and ", basename(svg_path), " -> ", out_dir)
}







# ------------------ READ OBJECT ------------------
message("Reading Seurat object: ", rds_path)
obj <- readRDS(rds_path)

# ------------------ EXPERIMENT DESIGN ------------------
untreated_lanes <- c(1, 3, 5, 7, 33, 35, 37, 39)
treated_lanes   <- c(2, 4, 6, 8, 34, 36, 38, 40)
lane_groups <- list(
  "untrained" = c(1, 2, 33, 34),
  "SYKi"      = c(3, 4, 35, 36),
  "MDP"       = c(5, 6, 37, 38),
  "BG"        = c(7, 8, 39, 40)
)
training_levels <- c(
  "untrained_Untreated","SYKi_Untreated","MDP_Untreated","BG_Untreated",
  "untrained_Treated","SYKi_Treated","MDP_Treated","BG_Treated"
)

stopifnot("Lane" %in% colnames(obj@meta.data))
stopifnot("predicted.celltype.l1" %in% colnames(obj@meta.data))

obj$Treatment <- ifelse(obj$Lane %in% untreated_lanes, "Untreated",
                        ifelse(obj$Lane %in% treated_lanes, "Treated", NA))
obj$LaneGroup <- dplyr::case_when(
  obj$Lane %in% lane_groups$untrained ~ "untrained",
  obj$Lane %in% lane_groups$SYKi      ~ "SYKi",
  obj$Lane %in% lane_groups$MDP       ~ "MDP",
  obj$Lane %in% lane_groups$BG        ~ "BG",
  TRUE ~ NA_character_
)
obj$LaneGroup_Treatment <- factor(paste(obj$LaneGroup, obj$Treatment, sep = "_"),
                                  levels = training_levels)

# ------------------ REDUCTION HANDLING ------------------
ensure_reduction <- function(obj) {
  if ("cca_clusters" %in% names(obj@reductions)) return(obj)
  for (nm in c("umap", "tsne")) {
    if (nm %in% names(obj@reductions)) {
      obj@reductions[["cca_clusters"]] <- obj@reductions[[nm]]
      message(sprintf("Using existing reduction '%s' as 'cca_clusters'", nm))
      return(obj)
    }
  }
  parent <- if ("integrated.cca" %in% names(obj@reductions)) "integrated.cca" else if ("pca" %in% names(obj@reductions)) "pca" else NA
  if (is.na(parent)) stop("No 2D reduction and no parent reduction ('integrated.cca' or 'pca') to compute from.")
  message("Computing UMAP from '", parent, "' → 'cca_clusters' ...")
  obj <- RunUMAP(obj, reduction = parent, dims = 1:20, reduction.name = "cca_clusters", verbose = FALSE)
  obj
}
obj <- ensure_reduction(obj)
emb_names <- colnames(obj@reductions$cca_clusters@cell.embeddings)
stopifnot(length(emb_names) >= 2)
dim1_name <- emb_names[1]
dim2_name <- emb_names[2]

# ------------------ LABEL CENTERS (for ggrepel) ------------------
centers <- obj@reductions$cca_clusters@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(".cell") %>%
  mutate(cell_type = obj$predicted.celltype.l1[.cell]) %>%
  group_by(cell_type) %>%
  summarise(
    !!dim1_name := median(.data[[dim1_name]], na.rm = TRUE),
    !!dim2_name := median(.data[[dim2_name]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(cell_type),
         nzchar(as.character(cell_type)),
         is.finite(.data[[dim1_name]]),
         is.finite(.data[[dim2_name]])) %>%
  distinct()

if (nrow(centers) > MAX_LABELS) {
  med1 <- median(centers[[dim1_name]]); med2 <- median(centers[[dim2_name]])
  centers <- centers %>%
    mutate(dist0 = ( .data[[dim1_name]] - med1 )^2 + ( .data[[dim2_name]] - med2 )^2) %>%
    arrange(dist0) %>% slice_head(n = MAX_LABELS) %>% select(-dist0)
}

# ------------------ MODULE SELECTION ------------------
all_module_cols <- grep("^ModuleScore", colnames(obj@meta.data), value = TRUE)
module_queries  <- unique(c(parse_cli_list(modules_arg), parse_file_list(modules_file)))
module_cols     <- if (length(module_queries)) pick_with_regex(all_module_cols, module_queries) else all_module_cols

if (length(module_queries) > 0) {
  if (length(module_cols) == 0) stop("No ModuleScore columns matched any of: ", paste(module_queries, collapse = ", "))
  message("Limiting to ", length(module_cols), " requested module(s).")
} else {
  message("No module filter provided; plotting all ", length(all_module_cols), " module(s).")
}

# ------------------ GENE SELECTION ------------------
genes_queries <- unique(c(parse_cli_list(genes_arg), parse_file_list(genes_file)))
all_genes_rna <- tryCatch(rownames(GetAssayData(obj, assay = "RNA")), error = function(e) character(0))

genes_to_plot <- character(0)
if (length(all_genes_rna)) {
  if (length(genes_queries)) {
    # accept exact matches and case-insensitive regex hits
    exact_g <- intersect(all_genes_rna, genes_queries)
    left_g  <- setdiff(genes_queries, exact_g)
    regex_g <- if (length(left_g)) grep(paste(left_g, collapse="|"), all_genes_rna, ignore.case = TRUE, value = TRUE) else character(0)
    genes_to_plot <- unique(c(exact_g, regex_g))
    if (length(genes_to_plot) == 0) {
      stop("No genes matched in RNA assay among: ", paste(genes_queries, collapse = ", "))
    } else {
      message("Plotting ", length(genes_to_plot), " requested gene(s).")
    }
  } else {
    # no gene filter provided; try defaults if present
    genes_to_plot <- intersect(DEFAULT_GENES, all_genes_rna)
    if (length(genes_to_plot)) {
      message("No gene list provided; defaulting to: ", paste(genes_to_plot, collapse = ", "))
    } else {
      message("No gene list provided and none of DEFAULT_GENES found in RNA; skipping gene plots.")
    }
  }
} else {
  message("RNA assay not available or has no rownames; skipping gene plots.")
}

# ------------------ COMMON PLOT BUILDERS ------------------
build_condition_panels <- function(feature_name, is_gene = FALSE, score_min = NULL, score_max = NULL) {
  # feature_name: if is_gene = TRUE, pass the RNA feature (e.g., "rna_GZMB"); else pass meta column name
  panels <- list()
  for (cond in training_levels) {
    cells_in_cond <- tryCatch(WhichCells(obj, expression = LaneGroup_Treatment == cond),
                              error = function(e) character(0))
    if (length(cells_in_cond) == 0) next

    p <- tryCatch({
      FeaturePlot(
        obj,
        features   = feature_name,
        reduction  = "cca_clusters",
        cells      = cells_in_cond,
        pt.size    = 0.5,
        min.cutoff = if (!is.null(score_min)) score_min else NA,
        max.cutoff = if (!is.null(score_max)) score_max else NA,
        cols       = c("lightblue", "blue", "red")
      ) +
        labs(title = NULL, subtitle = cond) +
        theme_minimal(base_size = 14)
    }, error = function(e) {
      message("FeaturePlot failed for ", cond, ": ", conditionMessage(e))
      NULL
    })

    if (!is.null(p) && ADD_LABELS && cond == "untrained_Untreated" && nrow(centers) > 0) {
      p <- tryCatch({
        p + ggrepel::geom_text_repel(
          data = centers,
          aes(x = .data[[dim1_name]], y = .data[[dim2_name]], label = cell_type),
          inherit.aes = FALSE,
          size = 3, fontface = "bold", color = "black",
          max.overlaps = 200, box.padding = 0.2, point.padding = 0.1,
          force = 0.5, min.segment.length = 0
        )
      }, error = function(e) {
        message("ggrepel skipped on labels: ", conditionMessage(e))
        p
      })
    }
    if (!is.null(p)) panels[[cond]] <- p
  }
  panels
}

build_violin <- function(y_col, title_y, dashed_uU = NA_real_, dashed_uT = NA_real_, show_text = TRUE) {
  v <- ggplot(obj@meta.data, aes(x = LaneGroup_Treatment, y = .data[[y_col]], fill = LaneGroup_Treatment)) +
    geom_violin(alpha = 0.5, trim = TRUE) +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
    scale_fill_manual(values = FILL_VALS) +
    labs(title = NULL, y = title_y, x = "Condition") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (is.finite(dashed_uU)) v <- v + annotate("segment", x = 1, xend = 4, y = dashed_uU, yend = dashed_uU,
                                              linetype = "dashed", color = "blue")
  if (is.finite(dashed_uT)) v <- v + annotate("segment", x = 5, xend = 8, y = dashed_uT, yend = dashed_uT,
                                              linetype = "dashed", color = "red")

  if (show_text) {
    med_by_group <- obj@meta.data %>%
      filter(!is.na(LaneGroup_Treatment)) %>%
      group_by(LaneGroup_Treatment) %>%
      summarise(med = median(.data[[y_col]], na.rm = TRUE), .groups = "drop") %>%
      mutate(label = format(round(med, MEDIAN_DIGITS), nsmall = MEDIAN_DIGITS))
    if (nrow(med_by_group) > 0) {
      y_const <- 1.25
      y_dat_min <- suppressWarnings(min(obj@meta.data[[y_col]], na.rm = TRUE))
      y_dat_max <- suppressWarnings(max(obj@meta.data[[y_col]], na.rm = TRUE))
      y_top     <- max(y_const, y_dat_max)
      pad       <- 0.05 * max(1, abs(y_top - y_dat_min))
      v <- v +
        geom_text(data = med_by_group,
                  aes(x = LaneGroup_Treatment, y = y_const, label = label),
                  inherit.aes = FALSE, size = 3.2, color = "black",
                  fontface = "bold", angle = 90, vjust = 0.5) +
        coord_cartesian(ylim = c(y_dat_min, y_top + pad), clip = "off")
    }
  }
  v
}

# ------------------ PLOT MODULES ------------------
if (length(module_cols)) {
  for (feature in module_cols) {
    if (!feature %in% colnames(obj@meta.data)) {
      message("Skipping missing module: ", feature)
      next
    }
    # scale for panels
    score_min <- min(obj@meta.data[[feature]], na.rm = TRUE)
    score_max <- stats::quantile(obj@meta.data[[feature]], 0.95, na.rm = TRUE)

    panels <- build_condition_panels(feature_name = feature, is_gene = FALSE,
                                     score_min = score_min, score_max = score_max)
    if (length(panels) == 0) {
      message("No condition panels built for module ", feature, " — skipping.")
      next
    }

    # dashed references from untrained medians
    untrained_meds <- obj@meta.data %>%
      filter(LaneGroup_Treatment %in% c("untrained_Untreated", "untrained_Treated")) %>%
      group_by(LaneGroup_Treatment) %>%
      summarise(median_score = median(.data[[feature]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = LaneGroup_Treatment, values_from = median_score)
    y_uU <- if ("untrained_Untreated" %in% names(untrained_meds)) untrained_meds$untrained_Untreated else NA_real_
    y_uT <- if ("untrained_Treated"   %in% names(untrained_meds)) untrained_meds$untrained_Treated   else NA_real_

    v <- build_violin(feature, "Module Score", dashed_uU = y_uU, dashed_uT = y_uT, show_text = SHOW_MEDIAN_TEXT)

    main_title <- paste("Combined Feature and Violin Plot:", feature)
    combined <- cowplot::plot_grid(
      cowplot::ggdraw() + cowplot::draw_label(main_title, fontface = 'bold', size = 16, hjust = 0, x = 0.02),
      cowplot::plot_grid(cowplot::plot_grid(plotlist = panels, ncol = 4, align = "v"), v, ncol = 1, rel_heights = c(2, 1)),
      ncol = 1, rel_heights = c(0.1, 1)
    )
    base <- paste0("plot_", feature)
    save_pdf_svg(combined, base, width = 12, height = 12)
  }
}

# ------------------ PLOT GENES ------------------
if (length(genes_to_plot)) {
  # Ensure we reference RNA features via "rna_<SYMBOL>" (Seurat v5 style)
  for (g in genes_to_plot) {
    gene_feat <- paste0("rna_", g)

    # pull expression for cutoffs
    expr_vec <- tryCatch(FetchData(obj, vars = gene_feat, slot = "data")[, 1], error = function(e) NULL)
    if (is.null(expr_vec)) {
      message("Skipping gene (no data): ", g)
      next
    }
    score_min <- suppressWarnings(quantile(expr_vec, 0.01, na.rm = TRUE))
    score_max <- if (toupper(g) == "GZMB") {
      suppressWarnings(quantile(expr_vec, 0.95, na.rm = TRUE))
    } else {
      suppressWarnings(quantile(expr_vec, 0.99, na.rm = TRUE))
    }

    panels <- build_condition_panels(feature_name = gene_feat, is_gene = TRUE,
                                     score_min = score_min, score_max = score_max)
    if (length(panels) == 0) {
      message("No condition panels built for gene ", g, " — skipping.")
      next
    }

    # attach expression to meta for violin
    expr_col <- paste0(g, "_expr")
    obj@meta.data[[expr_col]] <- tryCatch(FetchData(obj, vars = gene_feat)[, 1], error = function(e) NA_real_)

    # dashed medians from untrained groups (optionally non-zero)
    base_df <- obj@meta.data
    if (GENE_VIOLIN_NONZERO_ONLY) base_df <- base_df %>% filter(.data[[expr_col]] > 0)
    untrained_meds <- base_df %>%
      filter(LaneGroup_Treatment %in% c("untrained_Untreated", "untrained_Treated")) %>%
      group_by(LaneGroup_Treatment) %>%
      summarise(median_score = median(.data[[expr_col]], na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = LaneGroup_Treatment, values_from = median_score)
    y_uU <- if ("untrained_Untreated" %in% names(untrained_meds)) untrained_meds$untrained_Untreated else NA_real_
    y_uT <- if ("untrained_Treated"   %in% names(untrained_meds)) untrained_meds$untrained_Treated   else NA_real_

    # violin (match your previous styling; non-zero subset optional)
    violin_df <- if (GENE_VIOLIN_NONZERO_ONLY) obj@meta.data %>% filter(.data[[expr_col]] > 0) else obj@meta.data
    v <- ggplot(violin_df, aes(x = LaneGroup_Treatment, y = .data[[expr_col]], fill = LaneGroup_Treatment)) +
      geom_violin(alpha = 0.5) +
      geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA) +
      scale_fill_manual(values = FILL_VALS) +
      labs(title = paste("Violin Plot of", g, "Expression", if (GENE_VIOLIN_NONZERO_ONLY) "(Non-Zero)" else ""),
           y = "Expression Level", x = "Condition") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      { if (is.finite(y_uU)) annotate("segment", x = 1, xend = 4, y = y_uU, yend = y_uU,
                                      linetype = "dashed", color = "blue") else NULL } +
      { if (is.finite(y_uT)) annotate("segment", x = 5, xend = 8, y = y_uT, yend = y_uT,
                                      linetype = "dashed", color = "red") else NULL }

    if (SHOW_MEDIAN_TEXT) {
      all_meds <- violin_df %>%
        filter(!is.na(LaneGroup_Treatment)) %>%
        group_by(LaneGroup_Treatment) %>%
        summarise(median_score = median(.data[[expr_col]], na.rm = TRUE), .groups = "drop") %>%
        mutate(label = format(round(median_score, MEDIAN_DIGITS), nsmall = MEDIAN_DIGITS))
      if (nrow(all_meds) > 0) {
        y_const <- 1.25
        y_dat_min <- suppressWarnings(min(violin_df[[expr_col]], na.rm = TRUE))
        y_dat_max <- suppressWarnings(max(violin_df[[expr_col]], na.rm = TRUE))
        y_top     <- max(y_const, y_dat_max)
        pad       <- 0.05 * max(1, abs(y_top - y_dat_min))
        v <- v +
          geom_text(data = all_meds,
                    aes(x = LaneGroup_Treatment, y = y_const, label = label),
                    inherit.aes = FALSE, size = 3.2, color = "black",
                    fontface = "bold", angle = 90, vjust = 0.5) +
          coord_cartesian(ylim = c(y_dat_min, y_top + pad), clip = "off")
      }
    }

    main_title <- paste("Combined Feature and Violin Plot:", g)
    combined <- cowplot::plot_grid(
      cowplot::ggdraw() + cowplot::draw_label(main_title, fontface = 'bold', size = 16, hjust = 0, x = 0.02),
      cowplot::plot_grid(cowplot::plot_grid(plotlist = panels, ncol = 4, align = "v"), v, ncol = 1, rel_heights = c(2, 1)),
      ncol = 1, rel_heights = c(0.1, 1)
    )
    base <- paste0("plot_gene_", g)
    save_pdf_svg(combined, base, width = 12, height = 12)
  }
}

message("Done. Plots in: ", OUTDIR)
