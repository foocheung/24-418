


# =========================================================
# MESSY MERGED SCRIPT (WORKING VERSION)
# Combines:
# - ZNF peak2gene analysis
# - Locus-based accessibility
# - Ensembl overlap
# - Local vs non-local comparison
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
  library(biomaRt)
  library(GenomicRanges)
})

# =========================================================
# INPUT
# =========================================================



input_file <- "PEAKv3/v3_significant_peaks_BG_CD14_Mono.csv"
out_dir <- "ZNF_summary_outputs_FINAL"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE)

cat("Total rows:", nrow(df), "\n")

# =========================================================
# BASIC ZNF FILTER
# =========================================================
znf_all <- df %>% filter(grepl("^ZNF", geneName))
cat("ZNF rows:", nrow(znf_all), "\n")

# =========================================================
# SIGNIFICANCE
# =========================================================
znf_all <- znf_all %>%
  mutate(
    is_significant =
      !is.na(Correlation) &
      abs(Correlation) > 0.2 &
      !is.na(FDR.y) &
      FDR.y < 0.05 &
      !is.na(Log2FC) &
      abs(Log2FC) > 0.5 &
      !is.na(FDR.x) &
      FDR.x < 0.05
  )

znf_sig <- znf_all %>% filter(is_significant)
znf_not_sig <- znf_all %>% filter(!is_significant)

cat("Significant:", nrow(znf_sig), "\n")

# =========================================================
# SUMMARY
# =========================================================
summary_counts <- data.frame(
  total = nrow(znf_all),
  sig = nrow(znf_sig),
  nonsig = nrow(znf_not_sig)
)

write.csv(summary_counts, file.path(out_dir, "summary_counts.csv"), row.names = FALSE)

# =========================================================
# GENE SUMMARY
# =========================================================
znf_gene_summary <- znf_all %>%
  group_by(geneName) %>%
  summarise(
    total = n(),
    sig = sum(is_significant),
    mean_log2fc = mean(Log2FC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sig))

write.csv(znf_gene_summary, file.path(out_dir, "gene_summary.csv"), row.names = FALSE)

top_genes <- znf_gene_summary %>% filter(sig > 0) %>% pull(geneName)

# =========================================================
# PLOT 1
# =========================================================
p1 <- ggplot(summary_counts %>% pivot_longer(everything()),
             aes(name, value)) +
  geom_col() +
  theme_bw()

ggsave(file.path(out_dir, "plot_counts.pdf"), p1)

# =========================================================
# DOT PLOT
# =========================================================
plot_df <- znf_sig %>%
  filter(geneName %in% top_genes) %>%
  mutate(neglog10 = -log10(FDR.y))

p2 <- ggplot(plot_df,
             aes(Log2FC, geneName,
                 color = neglog10,
                 size = abs(Correlation))) +
  geom_point(alpha = 0.8) +
  theme_bw()

ggsave(file.path(out_dir, "dotplot.pdf"), p2, width = 10, height = 6)

# =========================================================
# LOCUS (PEAK SPAN METHOD)
# =========================================================
znf_loci <- znf_all %>%
  group_by(geneName, seqnames) %>%
  summarise(
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

write.csv(znf_loci, file.path(out_dir, "loci_span.csv"), row.names = FALSE)

# =========================================================
# UNIQUE PEAKS
# =========================================================
znf_peaks <- znf_all %>%
  distinct(seqnames, start, end, .keep_all = TRUE) %>%
  mutate(
    peak_DA = abs(Log2FC) > 0.5 & FDR.x < 0.05,
    peak_open = Log2FC > 0.5 & FDR.x < 0.05,
    peak_closed = Log2FC < -0.5 & FDR.x < 0.05
  )

write.csv(znf_peaks, file.path(out_dir, "unique_peaks.csv"), row.names = FALSE)

# =========================================================
# ENSEMBL COORDINATES
# =========================================================
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

znf_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = unique(znf_all$geneName),
  mart = mart
)

znf_coords <- znf_coords %>%
  filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name))

write.csv(znf_coords, file.path(out_dir, "ensembl_coords.csv"), row.names = FALSE)

# =========================================================
# GRANGES OVERLAP
# =========================================================
znf_gr <- GRanges(
  seqnames = znf_coords$chromosome_name,
  ranges = IRanges(
    start = znf_coords$start_position - 2000,
    end = znf_coords$end_position + 2000
  )
)

peaks_gr <- GRanges(
  seqnames = znf_peaks$seqnames,
  ranges = IRanges(start = znf_peaks$start, end = znf_peaks$end)
)

hits <- findOverlaps(peaks_gr, znf_gr)

peaks_overlap <- znf_peaks[queryHits(hits), ]

write.csv(peaks_overlap, file.path(out_dir, "ensembl_overlap.csv"), row.names = FALSE)

# =========================================================
# LOCAL VS NON-LOCAL
# =========================================================
znf_sig_unique <- znf_sig %>%
  distinct(seqnames, start, end, .keep_all = TRUE) %>%
  mutate(id = paste(seqnames, start, end))

znf_local <- peaks_overlap %>%
  mutate(id = paste(seqnames, start, end))

annot <- znf_sig_unique %>%
  mutate(
    locus = ifelse(id %in% znf_local$id, "local", "non-local"),
    neglog = -log10(FDR.y)
  )

write.csv(annot, file.path(out_dir, "local_vs_nonlocal.csv"), row.names = FALSE)

# =========================================================
# FINAL PLOT
# =========================================================
# ORDER GENES PROPERLY (by mean Log2FC)
gene_order <- annot %>%
  group_by(geneName) %>%
  summarise(mean_Log2FC = mean(Log2FC, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_Log2FC) %>%
  pull(geneName)

annot <- annot %>%
  mutate(geneName = factor(geneName, levels = gene_order))

# FINAL PLOT (ORDERED)
p3 <- ggplot(annot,
             aes(Log2FC, geneName,
                 color = neglog,
                 shape = locus)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.9) +
  theme_bw() +
  labs(
    title = "ZNF peaks: local vs non-local",
    x = "Log2FC",
    y = "ZNF gene"
  )

ggsave(file.path(out_dir, "local_vs_nonlocal_dot.pdf"), p3, width = 10, height = 8)
# =========================================================
# DONE
# =========================================================
cat("DONE. Output:", out_dir, "\n")
