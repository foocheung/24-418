#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(circlize)
  library(grid)
  library(ComplexHeatmap)  # Legend()
})

# -------------------------
# Inputs (your bin files)
# -------------------------
bins_cd4  <- "CD4_TCM_output_20250108_104840_FC/genome_bin_data_all_genes.csv"
bins_mono <- "CD14_mono_output_20241220_131510_FC/genome_bin_data_all_genes.csv"

# -------------------------
# Output directory
# -------------------------
out_dir <- "MS_V1/CIRCOS"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Helpers
# -------------------------
to_chr_prefix <- function(x) {
  x <- as.character(x)
  x <- gsub("^([0-9XYM]+)$", "chr\\1", x)
  sub("^chrMT$", "chrM", x)
}

read_bins <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)
  
  # Normalize first 3 columns to chr/start/end in case of unnamed extras
  stopifnot(ncol(df) >= 6)
  colnames(df)[1:3] <- c("chr", "start", "end")
  
  # Keep only needed columns if present
  keep <- c("chr", "start", "end", "Count_SKYi", "Count_MDP", "Count_BG")
  df <- df[, intersect(keep, names(df))]
  
  # Ensure all 3 datasets exist (error if not)
  missing_cols <- setdiff(c("Count_SKYi", "Count_MDP", "Count_BG"), names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }
  
  df %>%
    mutate(chr = to_chr_prefix(chr))
}

chrom_order_from <- function(df) {
  preferred <- paste0("chr", c(1:22, "X", "Y", "M"))
  present   <- unique(df$chr)
  ord <- preferred[preferred %in% present]
  if (length(ord) == 0) ord <- sort(present)
  ord
}

draw_bins_with_legend <- function(data, title_text = NULL) {
  dataset_colors <- c("Count_BG" = "red",
                      "Count_SKYi" = "blue",
                      "Count_MDP" = "green")
  
  chrom_order <- chrom_order_from(data)
  
  circos.clear()
  circos.par(cell.padding = c(0.01, 0.01, 0.01, 0.01))
  circos.initializeWithIdeogram(chromosome.index = chrom_order)
  
  for (dataset in names(dataset_colors)) {
    data_single <- data[, c("chr", "start", "end", dataset)]
    colnames(data_single) <- c("chr", "start", "end", "Count")
    data_single$chr <- factor(data_single$chr, levels = chrom_order)
    
    y_max <- max(data_single$Count, na.rm = TRUE)
    y_at  <- pretty(c(0, y_max), n = 5)
    first_sector <- TRUE
    
    circos.genomicTrack(
      data_single,
      ylim = c(0, y_max),
      panel.fun = function(region, value, ...) {
        circos.genomicRect(
          region  = region,
          value   = value,
          ybottom = 0,
          ytop    = value$Count,
          col     = dataset_colors[dataset],
          border  = dataset_colors[dataset],
          lwd     = 1
        )
        if (first_sector) {
          circos.yaxis(side = "left", at = y_at, labels = y_at, labels.cex = 0.5)
          first_sector <<- FALSE
        }
      }
    )
  }
  
  # Title + Legend
  if (!is.null(title_text)) {
    grid.text(title_text, x = unit(0.5, "npc"), y = unit(0.98, "npc"))
  }
  lgd <- Legend(
    labels    = c("BG", "SKYi", "MDP"),
    legend_gp = gpar(fill = c("red", "blue", "green"), col = c("red", "blue", "green")),
    type = "points", pch = 15, title = "Datasets"
  )
  draw(lgd, x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left","top"))
  
  circos.clear()
}

write_plots <- function(data, prefix, title_text) {
  pdf(file.path(out_dir, paste0(prefix, ".pdf")), width = 10, height = 10, useDingbats = FALSE)
  draw_bins_with_legend(data, title_text)
  dev.off()
  
  if (!requireNamespace("svglite", quietly = TRUE)) {
    install.packages("svglite", repos = "https://cloud.r-project.org")
  }
  svglite::svglite(file.path(out_dir, paste0(prefix, ".svg")), width = 10, height = 10)
  draw_bins_with_legend(data, title_text)
  dev.off()
}

# -------------------------
# Run for CD4 TCM bins
# -------------------------
data_cd4  <- read_bins(bins_cd4)
write_plots(data_cd4,  "CD4_TCM_FC_circos",  "CD4 TCM – FC bins")

# -------------------------
# Run for CD14 Mono bins
# -------------------------
data_mono <- read_bins(bins_mono)
write_plots(data_mono, "CD14_Mono_FC_circos", "CD14 Mono – FC bins")

message("Wrote circos plots to: ", normalizePath(out_dir))



################################




#!/usr/bin/env Rscript

# ---- Libraries ----
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(gridExtra)
library(UpSetR)
library(circlize)
library(IRanges)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg38)

library(GenomicRanges)  # (already loaded above)
library(readr)
library(circlize)

# ---- Data prep (unchanged) ----
mm2_hits <- read_csv("./mmc2.csv")  # Your insertion data
mm2_gr <- GRanges(
  seqnames = mm2_hits$Chromosome,
  ranges   = IRanges(start = mm2_hits$start, width = 1),
  strand   = mm2_hits$strand
)

mm2_gr_df <- as.data.frame(mm2_gr)

# (first version)
circos_data <- data.frame(
  chr   = as.character(mm2_gr_df$seqnames),
  start = mm2_gr_df$start,
  end   = mm2_gr_df$end
)

# (overwritten with midpoint version, as in your script)
mm2_gr_df <- as.data.frame(mm2_gr)
circos_data <- data.frame(
  chr      = as.character(mm2_gr_df$seqnames),
  position = (mm2_gr_df$start + mm2_gr_df$end) / 2
)

a <- read_csv("./CD4_PEAKS/t_significant_peaks_BG_CD4_TCM.csv") %>%
  dplyr::filter(!is.na(start)) %>%
  # dplyr::filter(FDR.x <= 0.5) %>%
  dplyr::filter(Log2FC > 1.5) %>%
  dplyr::filter(Correlation > 0.3) %>%
  dplyr::filter(FDR.y <= 0.05) %>%
  dplyr::select(chr = seqnames, start = start, end = end, value = Log2FC, FDR = FDR.y) %>%
  unique()

#data <- read_csv("./CD4_TCM_output_20241220_125350_FC/genome_bin_data_all_genes.csv")
data <- read_csv("./CD4_TCM_output_20250108_104840_FC/genome_bin_data_all_genes.csv")
colnames(data) <- c("chr", "start", "end", "Count_SKYi", "Count_MDP", "Count_BG")

# Only plotting Count_BG, per your script
dataset_colors <- c("Count_BG" = "green")

# Legend colors/labels (as in your script)
legend_colors <- c("HIV_Insertions" = "red", "Gene_Count_BG" = "green", "ATAC_Peaks_BG" = "blue")

# ---- Drawing routine (your plotting code, bundled) ----
draw_mmc_circos <- function() {
  circos.clear()
  
  # Initialize circos with ideogram
  circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:22, "X", "Y")))
  
  # Points track (insertion midpoints with random y)
  circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      circos.points(
        x = circos_data$position[circos_data$chr == CELL_META$sector.index],
        y = runif(sum(circos_data$chr == CELL_META$sector.index), min = 0, max = 1),
        pch = 16, cex = 0.5, col = "red"
      )
    }
  )
  
  # Binned BG counts as bars
  for (dataset in names(dataset_colors)) {
    data_single <- data[, c("chr", "start", "end", dataset)]
    colnames(data_single) <- c("chr", "start", "end", "Count")
    
    circos.genomicTrack(
      data_single,
      ylim = c(0, max(data_single$Count, na.rm = TRUE)),
      panel.fun = function(region, value, ...) {
        circos.genomicRect(
          region = region,
          value  = value,
          ybottom = 0,
          ytop    = value$Count,
          col     = dataset_colors[dataset],
          border  = dataset_colors[dataset],
          lwd     = 1
        )
      }
    )
  }
  
  # Legend (base graphics, as in your script)
  legend("topright", legend = names(legend_colors), fill = unname(legend_colors), title = "Datasets")
  
  # ATAC peaks (blue vertical lines)
  circos.genomicTrack(
    a,
    ylim = c(0, max(a$value, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value$value, col = "blue", type = "h", lwd = 2)
    }
  )
}

# ---- Render to PDF ----
pdf("MS_V1/CIRCOS/mmc_circos.pdf", width = 10, height = 10, useDingbats = FALSE)
draw_mmc_circos()
dev.off()

# ---- Render to SVG ----
if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite", repos = "https://cloud.r-project.org")
}
svglite::svglite("MS_V1/CIRCOS/mmc_circos.svg", width = 10, height = 10)
draw_mmc_circos()
dev.off()

