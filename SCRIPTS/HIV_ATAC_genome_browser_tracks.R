#!/usr/bin/env Rscript

############################################################
# HIV Genome Browser-Style ATAC-seq Track Visualization
#
# Purpose:
# This script visualizes HIV genomic annotations together with
# ATAC-seq accessibility signal across training conditions.
#
# Inputs:
# - sequence.fasta: HIV reference genome
# - sequence.gff3: HIV gene annotations
# - output.bed: ATAC-seq peak coordinates and scores
# - HIV_genes.info: optional HIV gene/scRNA-seq summary table
#
# Output:
# - Genome browser-style track plot saved as:
#     MS_V1/Genome_Browser/genome_browser_6.pdf
#     MS_V1/Genome_Browser/genome_browser_6.svg
#
# Notes:
# - This script is for visualization only.
# - It does not perform statistical testing.
############################################################


############################################################
# Libraries
############################################################

library(ggbio)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(ggrepel)
library(svglite)


############################################################
# Output Directory
############################################################

out_dir <- file.path("MS_V1", "Genome_Browser")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


############################################################
# Load HIV Genome, ATAC Peaks, and Gene Annotation
############################################################

hiv_genome <- import("sequence.fasta", format = "fasta")

peak_data <- read.table("output.bed", header = TRUE)

gff_file <- import("./sequence.gff3", format = "gff")
gene_annotations <- subset(gff_file, type == "gene")


############################################################
# Format ATAC Peak Table
############################################################

colnames(peak_data) <- c(
  "chromosome",
  "start",
  "end",
  "name",
  "score",
  "PMA_Activation",
  "Training_condition"
)

peak_data$start <- as.numeric(peak_data$start)
peak_data$end <- as.numeric(peak_data$end)
peak_data$score <- as.numeric(peak_data$score)


############################################################
# Convert ATAC Peaks to GRanges
############################################################

gr_peak_data <- GRanges(
  seqnames = Rle(peak_data$chromosome),
  ranges = IRanges(
    start = peak_data$start,
    end = peak_data$end
  ),
  PMA_Activation = peak_data$PMA_Activation,
  Training_condition = peak_data$Training_condition,
  score = peak_data$score
)


############################################################
# Helper Function: Summarize ATAC Signal Across Genome
#
# For each genomic position, this function sums the scores of
# all overlapping ATAC peaks. The resulting table is used to
# create genome browser-style accessibility tracks.
############################################################

create_cumulative_scores <- function(subset_atac) {
  
  df <- as.data.frame(subset_atac)
  
  total_counts <- data.frame(
    start = integer(),
    total_count = numeric()
  )
  
  for (pos in seq(min(df$start), max(df$end))) {
    overlapping_rows <- df[df$start <= pos & df$end > pos, ]
    total_count <- sum(overlapping_rows$score)
    
    total_counts <- rbind(
      total_counts,
      data.frame(
        start = pos,
        total_count = total_count
      )
    )
  }
  
  total_counts$cumulative_count <- cumsum(total_counts$total_count)
  
  total_counts
}


############################################################
# Build ATAC-seq Tracks by Condition and Treatment
############################################################

atac_tracks <- list()
max_y_limit <- 0

for (group in unique(peak_data$Training_condition)) {
  for (treatment in unique(peak_data$PMA_Activation)) {
    
    subset_atac <- gr_peak_data[
      gr_peak_data$Training_condition == group &
        gr_peak_data$PMA_Activation == treatment
    ]
    
    if (length(subset_atac) == 0) next
    
    cumulative_data <- create_cumulative_scores(subset_atac)
    
    max_y_limit <- max(
      max_y_limit,
      max(cumulative_data$total_count, na.rm = TRUE)
    )
    
    track_name <- paste(group, treatment)
    
    atac_tracks[[track_name]] <- ggplot(cumulative_data) +
      geom_col(
        aes(x = start, y = total_count),
        fill = "black",
        alpha = 0.7
      ) +
      labs(
        title = paste("ATAC-seq:", group, "-", treatment),
        x = "Genomic Position",
        y = NULL
      ) +
      theme_minimal() +
      xlim(0, max(cumulative_data$start) + 1000) +
      ylim(0, max_y_limit + 1) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid.major.y = element_line(color = "gray", linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_line(linewidth = 0.5)
      )
  }
}


############################################################
# Optional scRNA-seq HIV Gene Track Preparation
#
# This section reads HIV gene count information and converts
# it into GRanges format. In the final manuscript plot below,
# scRNA-seq tracks are currently not included, but this object
# is retained for optional future plotting.
############################################################

scrna_data <- read.table("./HIV_genes.info", header = TRUE, sep = "\t")

gene_positions <- data.frame(
  gene = c("HIV1gp1", "HIV1gp3", "HIV1gp5", "HIV1gp8"),
  start = c(336, 4587, 5377, 5771),
  end = c(4642, 5165, 7970, 8341)
)

merged_data <- merge(
  scrna_data,
  gene_positions,
  by.x = "Gene",
  by.y = "gene"
)

scRNAseq_gr <- GRanges(
  seqnames = Rle("NC_001802.1"),
  ranges = IRanges(
    start = merged_data$start,
    end = merged_data$end
  ),
  gene = merged_data$Gene,
  type = "gene",
  expression = merged_data$Count,
  group = merged_data$Group,
  treatment = merged_data$Treatment,
  genome = merged_data$Genome
)


############################################################
# Helper Function: Add Track if Present
#
# If a requested track exists, it is returned.
# If not, an empty placeholder track is returned so the final
# plot layout does not fail.
############################################################

add_track_if_exists <- function(track_list, track_name, xlim_range = NULL) {
  
  if (track_name %in% names(track_list)) {
    
    track <- track_list[[track_name]] +
      theme(legend.position = "none")
    
    if (!is.null(xlim_range)) {
      track <- track + xlim(xlim_range[1], xlim_range[2])
    }
    
    return(track)
  }
  
  empty_track <- ggplot() +
    geom_blank() +
    theme_void() +
    labs(title = track_name)
  
  if (!is.null(xlim_range)) {
    empty_track <- empty_track + xlim(xlim_range[1], xlim_range[2])
  }
  
  empty_track
}


############################################################
# Build Optional scRNA-seq Tracks
#
# These tracks are created but not currently included in the
# final plot. They can be added back into the tracks() call
# if scRNA-seq gene-level visualization is needed.
############################################################

scRNAseq_tracks <- list()

max_stack <- max(
  table(scRNAseq_gr$group, scRNAseq_gr$treatment),
  na.rm = TRUE
)

for (group in unique(merged_data$Group)) {
  for (treatment in unique(merged_data$Treatment)) {
    
    subset_scrna <- scRNAseq_gr[
      scRNAseq_gr$group == group &
        scRNAseq_gr$treatment == treatment
    ]
    
    if (length(subset_scrna) == 0) next
    
    df <- as.data.frame(subset_scrna)
    df <- df[order(df$start), ]
    
    df$stack_y <- ave(
      df$start,
      df$seqnames,
      FUN = function(x) seq_along(x)
    )
    
    track_name <- paste(group, treatment)
    
    scRNAseq_tracks[[track_name]] <- ggplot(df) +
      geom_tile(
        aes(
          x = (start + end) / 2,
          y = stack_y,
          width = end - start,
          height = 0.8
        ),
        color = "black"
      ) +
      labs(
        title = paste("scRNA-seq:", group, "-", treatment),
        y = NULL
      ) +
      scale_y_continuous(labels = NULL) +
      theme_minimal() +
      ylim(0, max_stack + 1) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_blank()
      )
  }
}


############################################################
# HIV Gene Annotation Track
############################################################

xlim_range <- c(
  min(start(gene_annotations)),
  max(end(gene_annotations))
)

gene_annotations$gene_id <- gene_annotations$Name

gene_annotation_track <- autoplot(gene_annotations) +
  geom_text_repel(
    aes(
      x = (start(gene_annotations) + end(gene_annotations)) / 2,
      y = -0.5,
      label = gene_annotations$gene_id
    ),
    size = 5,
    fontface = "bold",
    nudge_x = 0.005,
    direction = "x",
    segment.color = NA,
    max.overlaps = 5
  ) +
  labs(title = "HIV Genes with IDs") +
  theme_minimal() +
  scale_x_continuous(limits = c(xlim_range[1], xlim_range[2]))


############################################################
# Assemble Final Genome Browser Plot
#
# The final figure includes:
# - HIV gene annotation track
# - PMA-treated ATAC-seq tracks for:
#     untrained, SYKi, MDP, BG
#
# Untreated and scRNA-seq tracks are available above but are
# not included in this final version.
############################################################

p <- tracks(
  Genes = gene_annotation_track + theme(legend.position = "none"),
  "ATAC\nUntrained\nTreated" = add_track_if_exists(
    atac_tracks,
    "untrained Treated",
    xlim_range
  ),
  "ATAC\nSYKi\nTreated" = add_track_if_exists(
    atac_tracks,
    "SYKi Treated",
    xlim_range
  ),
  "ATAC\nMDP\nTreated" = add_track_if_exists(
    atac_tracks,
    "MDP Treated",
    xlim_range
  ),
  "ATAC\nBG\nTreated" = add_track_if_exists(
    atac_tracks,
    "BG Treated",
    xlim_range
  ),
  title = "",
  label.width = unit(4.5, "lines")
)


############################################################
# Save Final Plot
############################################################

final_plot <- p +
  theme_tracks_sunset() +
  guides(fill = "none")

pdf_path <- file.path(out_dir, "genome_browser_6.pdf")
svg_path <- file.path(out_dir, "genome_browser_6.svg")

grDevices::pdf(pdf_path, width = 12, height = 5)
print(final_plot)
dev.off()

svglite::svglite(svg_path, width = 12, height = 5)
print(final_plot)
dev.off()

message("Saved genome browser plots to: ", out_dir)
