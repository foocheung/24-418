# Load necessary libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(gridExtra)
library(circlize)
library(ggplot2)
library(grid)
library(gtable)
library(svglite)  # NEW

# === OUTPUT FOLDERS ===
# Create canonical output folder for this figure set
out_dir <- file.path("MS_V1", "ATAC_pathways")  # NEW
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)  # NEW

# (optional: keep timestamped folder if you still want it)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- paste0("output_", timestamp, "_FDR_FC")
dir.create(output_dir, showWarnings = FALSE)

# Function to process gene data and perform analysis
process_genes <- function(file_path, output_file_prefix, plot_title_prefix) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  data <- read_csv(file_path)
  
  # Filter correlated and anti-correlated genes
  correlated_genes <- data %>%
    #filter(FDR.x <= 0.05, !is.na(start), Log2FC > 1.5, Correlation > 0.3, FDR.y <= 0.05) %>%
    filter(!is.na(start), Log2FC > 1.5, Correlation > 0.3, FDR.y <= 0.05) %>%
    distinct(geneName) %>% pull(geneName)
  
  anti_correlated_genes <- data %>%
    filter(!is.na(start), Log2FC > 1.5, Correlation < -0.3, FDR.y <= 0.05) %>%
    # filter(FDR.x <= 0.05, !is.na(start), Log2FC > 1.5, Correlation < -0.3, FDR.y <= 0.05) %>%
    distinct(geneName) %>% pull(geneName)
  
  # Helper function for pathway analysis
  analyze_and_plot <- function(genes, plot_title, tag) {
    if (length(genes) == 0) return(NULL)
    
    entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (is.null(entrez_ids) || nrow(entrez_ids) == 0) return(NULL)
    
    reactome_results <- enrichPathway(gene = entrez_ids$ENTREZID, organism = "human",
                                      pvalueCutoff = 0.05, readable = TRUE)
    go_results <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID", ont = "ALL",
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    
    reactome_plot <- if (nrow(as.data.frame(reactome_results)) > 0) {
      barplot(reactome_results, showCategory = 20, title = paste(plot_title, "Reactome"))
    } else { NULL }
    
    go_plot <- if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
      dotplot(go_results, showCategory = 20, title = paste(plot_title, "GO Enrichment")) +
        theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
              axis.text = element_text(size = 12),
              plot.title = element_text(hjust = 0.5))
    } else { NULL }
    
    # --- NEW: also write individual SVGs when plots exist ---
    if (!is.null(reactome_plot)) {
      ggsave(filename = file.path(out_dir, paste0(tag, "_Reactome.svg")),  # NEW
             plot = reactome_plot, width = 8, height = 6, units = "in")    # NEW
    }
    if (!is.null(go_plot)) {
      ggsave(filename = file.path(out_dir, paste0(tag, "_GO.svg")),        # NEW
             plot = go_plot, width = 8, height = 6, units = "in")          # NEW
    }
    # --------------------------------------------------------
    
    return(list(reactome_plot = reactome_plot, go_plot = go_plot))
  }
  
  correlated_plots <- analyze_and_plot(
    correlated_genes, paste0(plot_title_prefix, " Correlated"),
    tag = paste0(output_file_prefix, "_correlated"))  # NEW tag
  
  anti_correlated_plots <- analyze_and_plot(
    anti_correlated_genes, paste0(plot_title_prefix, " Anti-Correlated"),
    tag = paste0(output_file_prefix, "_anti"))  # NEW tag
  
  return(list(correlated = correlated_plots, anti_correlated = anti_correlated_plots))
}

# File paths and output configurations
file_paths <- c(
  "PEAKv3/v3_significant_peaks_BG_CD14_Mono.csv", 
  "PEAKv3/v3_significant_peaks_MDP_CD14_Mono.csv", 
  "PEAKv3/v3_significant_peaks_SYKi_CD14_Mono.csv"
)

output_file_prefixes <- sub("PEAKv3/v3_significant_peaks_", "", file_paths)
output_file_prefixes <- sub(".csv", "", output_file_prefixes)
plot_title_prefixes <- output_file_prefixes

# Process files and generate plots
all_plots <- lapply(seq_along(file_paths), function(i) {
  process_genes(file_paths[i], output_file_prefixes[i], plot_title_prefixes[i])
})

# PDF (unchanged except path)
pdf(file.path(out_dir, "4_pathway_figure.pdf"), height = 12, width = 12)  # NEW path
for (i in seq_along(all_plots)) {
  cat("\nSaving plots for dataset:", i, "\n")
  correlated_reactome <- all_plots[[i]]$correlated$reactome_plot
  correlated_go <- all_plots[[i]]$correlated$go_plot
  anti_correlated_reactome <- all_plots[[i]]$anti_correlated$reactome_plot
  anti_correlated_go <- all_plots[[i]]$anti_correlated$go_plot
  
  plot_list <- Filter(Negate(is.null), list(
    correlated_reactome, correlated_go
    # anti_correlated_reactome, anti_correlated_go
  ))
  
  if (length(plot_list) > 0) {
    grid.arrange(grobs = plot_list, ncol = 2,
                 top = paste("Pathway Analysis for Dataset", i))
    
    # --- NEW: save same arranged panel as SVG ---
    svglite::svglite(file.path(out_dir, paste0("4_pathway_figure_dataset", i, ".svg")),  # NEW
                     width = 12, height = 12)                                            # NEW
    grid.arrange(grobs = plot_list, ncol = 2,
                 top = paste("Pathway Analysis for Dataset", i))                          # NEW
    dev.off()                                                                             # NEW
    # -------------------------------------------
  } else {
    message("No valid plots for dataset ", i)
  }
}
dev.off()
message("All plots saved to PDF and per-dataset SVGs in MS_V1/ATAC_pathways")  # NEW

# === Big combined panel (correlated only) ===
base_font_size <- 16
correlated_plot_list <- list()

for (i in seq_along(all_plots)) {
  reactome_plot <- all_plots[[i]]$correlated$reactome_plot
  go_plot <- all_plots[[i]]$correlated$go_plot
  
  if (!is.null(reactome_plot)) {
    reactome_plot <- reactome_plot +
      ggtitle(paste("Dataset", i, "- Reactome")) +
      theme_minimal(base_size = base_font_size) +
      theme(plot.margin = margin(15, 10, 15, 20))
    g <- ggplotGrob(reactome_plot)
    g$layout$clip[g$layout$name == "panel"] <- "off"
    correlated_plot_list <- append(correlated_plot_list, list(g))
  }
  if (!is.null(go_plot)) {
    go_plot <- go_plot +
      ggtitle(paste("Dataset", i, "- GO")) +
      theme_minimal(base_size = base_font_size) +
      theme(plot.margin = margin(15, 10, 15, 20))
    g <- ggplotGrob(go_plot)
    g$layout$clip[g$layout$name == "panel"] <- "off"
    correlated_plot_list <- append(correlated_plot_list, list(g))
  }
}

if (length(correlated_plot_list) > 0) {
  # PDF (original)
  pdf(file.path(out_dir, "MS_FIGURE_mono_correlated_pathways_all_datasets.pdf"), height = 33, width = 13)  # NEW path
  do.call(grid.arrange, c(correlated_plot_list, ncol = 2,
                          top = "Correlated Pathway Plots Across All Datasets"))
  dev.off()
  
  # --- NEW: SVG of the same combined figure ---
  svglite::svglite(file.path(out_dir, "MS_FIGURE_mono_correlated_pathways_all_datasets.svg"),
                   width = 13, height = 33)
  do.call(grid.arrange, c(correlated_plot_list, ncol = 2,
                          top = "Correlated Pathway Plots Across All Datasets"))
  dev.off()
  
  message("✅ Saved combined correlated plots to PDF and SVG in MS_V1/ATAC_pathways")  # NEW
} else {
  message("⚠️ No correlated plots found.")
}



####################################################################

# # Load necessary libraries
# library(tidyverse)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
# library(gridExtra)
# library(circlize)
# library(ggplot2)
# library(grid)
# library(gtable)
# library(svglite)  # NEW
# 
# # Output folder (reuse same)
# out_dir <- file.path("MS_V1", "ATAC_pathways")        # NEW
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)  # NEW
# 
# # (optional) keep timestamped output_dir if useful elsewhere
# timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
# output_dir <- paste0("output_", timestamp, "_FDR_FC")
# dir.create(output_dir, showWarnings = FALSE)
# 
# process_genes <- function(file_path, output_file_prefix, plot_title_prefix) {
#   if (!file.exists(file_path)) stop(paste("File not found:", file_path))
#   data <- read_csv(file_path)
#   
#   correlated_genes <- data %>%
#     filter(!is.na(start), Log2FC > 1.5, Correlation > 0.3, FDR.y <= 0.05) %>%
#     distinct(geneName) %>% pull(geneName)
#   
#   anti_correlated_genes <- data %>%
#     filter(!is.na(start), Log2FC > 1.5, Correlation < -0.3, FDR.y <= 0.05) %>%
#     distinct(geneName) %>% pull(geneName)
#   
#   analyze_and_plot <- function(genes, plot_title, tag) {
#     if (length(genes) == 0) return(NULL)
#     entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#     if (is.null(entrez_ids) || nrow(entrez_ids) == 0) return(NULL)
#     
#     reactome_results <- enrichPathway(gene = entrez_ids$ENTREZID, organism = "human",
#                                       pvalueCutoff = 0.05, readable = TRUE)
#     go_results <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db,
#                            keyType = "ENTREZID", ont = "ALL",
#                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
#     
#     reactome_plot <- if (nrow(as.data.frame(reactome_results)) > 0) {
#       barplot(reactome_results, showCategory = 20, title = paste(plot_title, "Reactome"))
#     } else { NULL }
#     
#     go_plot <- if (!is.null(go_results) && nrow(as.data.frame(go_results)) > 0) {
#       dotplot(go_results, showCategory = 20, title = paste(plot_title, "GO Enrichment")) +
#         theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
#               axis.text = element_text(size = 12),
#               plot.title = element_text(hjust = 0.5))
#     } else { NULL }
#     
#     # NEW: write individual SVGs for each plot
#     if (!is.null(reactome_plot)) {
#       ggsave(filename = file.path(out_dir, paste0(tag, "_Reactome.svg")),
#              plot = reactome_plot, width = 8, height = 6, units = "in")
#     }
#     if (!is.null(go_plot)) {
#       ggsave(filename = file.path(out_dir, paste0(tag, "_GO.svg")),
#              plot = go_plot, width = 8, height = 6, units = "in")
#     }
#     
#     return(list(reactome_plot = reactome_plot, go_plot = go_plot))
#   }
#   
#   correlated_plots <- analyze_and_plot(
#     correlated_genes, paste0(plot_title_prefix, " Correlated"),
#     tag = paste0(output_file_prefix, "_correlated"))
#   anti_correlated_plots <- analyze_and_plot(
#     anti_correlated_genes, paste0(plot_title_prefix, " Anti-Correlated"),
#     tag = paste0(output_file_prefix, "_anti"))
#   
#   return(list(correlated = correlated_plots, anti_correlated = anti_correlated_plots))
# }
# 
# # File paths and output configurations
# file_paths <- c(
#   "PEAKv3/v3_significant_peaks_BG_CD4_TCM.csv", 
#   "PEAKv3/v3_significant_peaks_MDP_CD4_TCM.csv", 
#   "PEAKv3/v3_significant_peaks_SYKi_CD4_TCM.csv"
# )
# 
# output_file_prefixes <- sub("PEAKv3/v3_significant_peaks_", "", file_paths)
# output_file_prefixes <- sub(".csv", "", output_file_prefixes)
# plot_title_prefixes <- output_file_prefixes
# 
# all_plots <- lapply(seq_along(file_paths), function(i) {
#   process_genes(file_paths[i], output_file_prefixes[i], plot_title_prefixes[i])
# })

# PDF (unchanged except path)
# pdf(file.path(out_dir, "4_pathway_figure_CD4TCM.pdf"), height = 12, width = 12)  # NEW path
# for (i in seq_along(all_plots)) {
#   cat("\nSaving plots for dataset:", i, "\n")
#   correlated_reactome <- all_plots[[i]]$correlated$reactome_plot
#   correlated_go <- all_plots[[i]]$correlated$go_plot
#   anti_correlated_reactome <- all_plots[[i]]$anti_correlated$reactome_plot
#   anti_correlated_go <- all_plots[[i]]$anti_correlated$go_plot
#   
#   plot_list <- Filter(Negate(is.null), list(
#     correlated_reactome, correlated_go
#     # anti_correlated_reactome, anti_correlated_go
#   ))
#   
#   if (length(plot_list) > 0) {
#     grid.arrange(grobs = plot_list, ncol = 2,
#                  top = paste("Pathway Analysis for Dataset", i))
#     
#     # NEW: SVG per-dataset arranged panel
#     svglite::svglite(file.path(out_dir, paste0("4_pathway_figure_CD4TCM_dataset", i, ".svg")),
#                      width = 12, height = 12)
#     grid.arrange(grobs = plot_list, ncol = 2,
#                  top = paste("Pathway Analysis for Dataset", i))
#     dev.off()
#   } else {
#     message("No valid plots for dataset ", i)
#   }
# }
# dev.off()
# message("All CD4 TCM plots saved to PDF and per-dataset SVGs in MS_V1/ATAC_pathways")  # NEW

# === Big combined panel (correlated only) ===
base_font_size <- 16
correlated_plot_list <- list()

for (i in seq_along(all_plots)) {
  reactome_plot <- all_plots[[i]]$correlated$reactome_plot
  go_plot <- all_plots[[i]]$correlated$go_plot
  
  if (!is.null(reactome_plot)) {
    reactome_plot <- reactome_plot +
      ggtitle(paste("Dataset", i, "- Reactome")) +
      theme_minimal(base_size = base_font_size) +
      theme(plot.margin = margin(15, 10, 15, 20))
    g <- ggplotGrob(reactome_plot)
    g$layout$clip[g$layout$name == "panel"] <- "off"
    correlated_plot_list <- append(correlated_plot_list, list(g))
  }
  if (!is.null(go_plot)) {
    go_plot <- go_plot +
      ggtitle(paste("Dataset", i, "- GO")) +
      theme_minimal(base_size = base_font_size) +
      theme(plot.margin = margin(15, 10, 15, 20))
    g <- ggplotGrob(go_plot)
    g$layout$clip[g$layout$name == "panel"] <- "off"
    correlated_plot_list <- append(correlated_plot_list, list(g))
  }
}

if (length(correlated_plot_list) > 0) {
  # PDF (original)
  pdf(file.path(out_dir, "FIGURE_CD4T_correlated_pathways_all_datasets.pdf"),
      height = 35, width = 13)  # NEW path
  do.call(grid.arrange, c(correlated_plot_list, ncol = 2,
                          top = "Correlated Pathway Plots Across All Datasets"))
  dev.off()
  
  # NEW: SVG of the same combined figure
  svglite::svglite(file.path(out_dir, "FIGURE_CD4T_correlated_pathways_all_datasets.svg"),
                   width = 13, height = 35)
  do.call(grid.arrange, c(correlated_plot_list, ncol = 2,
                          top = "Correlated Pathway Plots Across All Datasets"))
  dev.off()
  
  message("✅ Saved CD4 TCM combined correlated plots to PDF and SVG in MS_V1/ATAC_pathways")  # NEW
} else {
  message("⚠️ No correlated plots found.")
}


