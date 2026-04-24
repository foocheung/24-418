############################################################
# ATAC-seq ZNF Peak2Gene vs Local Locus Analysis
#
# Purpose:
# This script distinguishes between:
# 1) Peaks linked to ZNF genes (peak2gene relationships)
# 2) Peaks physically located at ZNF loci (gene body ± 2 kb)
#
# Key steps:
# - Filter peak2gene table for ZNF-linked genes
# - Apply significance thresholds (correlation, FDR, Log2FC)
# - Deduplicate peaks while retaining strongest peak-gene link
# - Retrieve genomic coordinates for ZNF genes (Ensembl)
# - Identify peaks overlapping ZNF loci (±2 kb)
#
# Outputs:
# - Peak2gene ZNF summaries
# - Local vs non-local peak classification
# - Gene-level summaries
# - Visualization of filtering progression
# - BED and GFF exports for local peaks
#
# Interpretation:
# - Peak2gene links capture regulatory associations (can be distal)
# - Local peaks represent direct chromatin accessibility at ZNF loci
#
# Notes:
# - This script combines regulatory inference (peak2gene)
#   with physical genomic overlap (locus-based filtering)
# - Results should distinguish distal regulation vs local accessibility
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(biomaRt)
  library(GenomicRanges)
  library(patchwork)
  library(readr)
})

# =========================================================
# INPUTS
# =========================================================
# If df is already loaded in memory, this script will use it.
# Otherwise uncomment the read.csv line below.

peak2gene_file <- "./CD4_PEAKS/t_significant_peaks_BG_CD14_Mono.csv"

out_dir <- "ZNF_peak2gene_clear_filters"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Settings
znf_pattern      <- "^ZNF"
cor_cutoff       <- 0.2
link_fdr_cutoff  <- 0.05
peak_fdr_cutoff  <- 0.05
log2fc_cutoff    <- 0.5
flank_bp         <- 2000
top_n            <- 20

# =========================================================
# 1. LOAD PEAK2GENE TABLE
# =========================================================
# df <- read.csv(peak2gene_file, stringsAsFactors = FALSE)

if (!exists("df")) {
  stop("Object 'df' not found. Either load it first or uncomment read.csv().")
}

cat("Total rows in input peak2gene table:", nrow(df), "\n")

required_cols <- c("geneName", "seqnames", "start", "end", "Correlation", "FDR.y")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

has_peak_fdr <- "FDR.x" %in% colnames(df)
has_log2fc   <- "Log2FC" %in% colnames(df)

# =========================================================
# HELPER: RANK ROWS SO DISTINCT() KEEPS THE MOST SIGNIFICANT
# =========================================================
# For duplicate peak coordinates, rank rows by:
# 1) smallest link FDR
# 2) largest absolute correlation
# 3) largest absolute Log2FC if available
rank_peak_links <- function(dat) {
  dat2 <- dat %>%
    mutate(
      rank_FDRy = ifelse(is.na(FDR.y), Inf, FDR.y),
      rank_absCorr = ifelse(is.na(Correlation), -Inf, abs(Correlation))
    )
  
  if (has_log2fc) {
    dat2 <- dat2 %>%
      mutate(rank_absLog2FC = ifelse(is.na(Log2FC), -Inf, abs(Log2FC))) %>%
      arrange(rank_FDRy, desc(rank_absCorr), desc(rank_absLog2FC))
  } else {
    dat2 <- dat2 %>%
      arrange(rank_FDRy, desc(rank_absCorr))
  }
  
  dat2 %>%
    dplyr::select(-starts_with("rank_"))
}

# =========================================================
# 2. FILTER 1 = KEEP PEAK2GENE ROWS LINKED TO ZNF GENES
# =========================================================
# This does NOT mean the peak is physically near the ZNF gene.
# It only means ArchR linked the peak to a gene whose name starts with ZNF.

znf_p2g_all <- df %>%
  filter(grepl(znf_pattern, geneName))

cat("Rows linked to ZNF genes:", nrow(znf_p2g_all), "\n")

write.csv(
  znf_p2g_all,
  file.path(out_dir, "01_ZNF_peak2gene_all_rows.csv"),
  row.names = FALSE
)

# =========================================================
# 3. FILTER 2 = MARK WHICH ZNF-LINKED PEAK2GENE ROWS ARE SIGNIFICANT
# =========================================================
# This is still peak2gene space, not physical locus overlap yet.

znf_p2g_all <- znf_p2g_all %>%
  mutate(
    pass_corr     = !is.na(Correlation) & abs(Correlation) > cor_cutoff,
    pass_link_fdr = !is.na(FDR.y) & FDR.y < link_fdr_cutoff,
    pass_peak_fdr = if (has_peak_fdr) !is.na(FDR.x) & FDR.x < peak_fdr_cutoff else TRUE,
    pass_log2fc   = if (has_log2fc) !is.na(Log2FC) & abs(Log2FC) > log2fc_cutoff else TRUE,
    is_significant = pass_corr & pass_link_fdr & pass_peak_fdr & pass_log2fc
  )

znf_p2g_sig <- znf_p2g_all %>% filter(is_significant)
znf_p2g_nonsig <- znf_p2g_all %>% filter(!is_significant)

cat("Significant ZNF-linked peak2gene rows:", nrow(znf_p2g_sig), "\n")
cat("Non-significant ZNF-linked peak2gene rows:", nrow(znf_p2g_nonsig), "\n")

write.csv(
  znf_p2g_sig,
  file.path(out_dir, "02_ZNF_peak2gene_significant_rows.csv"),
  row.names = FALSE
)

write.csv(
  znf_p2g_nonsig,
  file.path(out_dir, "03_ZNF_peak2gene_nonsignificant_rows.csv"),
  row.names = FALSE
)

# =========================================================
# 4. SUMMARIZE THE PEAK2GENE ZNF VIEW
# =========================================================
znf_p2g_counts <- data.frame(
  step = c(
    "All ZNF-linked peak2gene rows",
    "Significant ZNF-linked peak2gene rows",
    "Non-significant ZNF-linked peak2gene rows"
  ),
  count = c(
    nrow(znf_p2g_all),
    nrow(znf_p2g_sig),
    nrow(znf_p2g_nonsig)
  )
)

write.csv(
  znf_p2g_counts,
  file.path(out_dir, "04_ZNF_peak2gene_counts.csv"),
  row.names = FALSE
)

znf_gene_summary <- znf_p2g_all %>%
  group_by(geneName) %>%
  summarise(
    total_rows = n(),
    sig_rows = sum(is_significant, na.rm = TRUE),
    nonsig_rows = sum(!is_significant, na.rm = TRUE),
    percent_sig = 100 * sig_rows / total_rows,
    mean_corr = mean(Correlation, na.rm = TRUE),
    min_link_fdr = min(FDR.y, na.rm = TRUE),
    mean_log2fc = if (has_log2fc) mean(Log2FC, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(desc(sig_rows), desc(percent_sig), min_link_fdr)

write.csv(
  znf_gene_summary,
  file.path(out_dir, "05_ZNF_peak2gene_gene_summary.csv"),
  row.names = FALSE
)

top_znf_genes <- znf_gene_summary %>%
  filter(sig_rows > 0) %>%
  slice_head(n = top_n) %>%
  pull(geneName)

# =========================================================
# 5. PLOT A = ZNF PEAK2GENE VIEW
# =========================================================
p_counts_before <- ggplot(znf_p2g_counts, aes(x = step, y = count, fill = step)) +
  geom_col(show.legend = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(
    title = "ZNF peak2gene view",
    subtitle = "Rows where the linked gene name starts with ZNF",
    x = "",
    y = "Number of rows"
  )

ggsave(
  file.path(out_dir, "Figure_01_ZNF_peak2gene_counts.pdf"),
  p_counts_before,
  width = 8,
  height = 5
)

plot_before_df <- znf_p2g_sig %>%
  filter(geneName %in% top_znf_genes) %>%
  mutate(
    geneName = factor(geneName, levels = rev(top_znf_genes)),
    neglog10_linkFDR = -log10(FDR.y)
  )

if (nrow(plot_before_df) > 0) {
  p_before_dot <- ggplot(
    plot_before_df,
    aes(
      x = if (has_log2fc) Log2FC else Correlation,
      y = geneName,
      size = abs(Correlation),
      color = neglog10_linkFDR
    )
  ) +
    geom_point(alpha = 0.85) +
    theme_bw(base_size = 12) +
    labs(
      title = "ZNF-linked significant peak2gene peaks",
      subtitle = "These peaks are linked to ZNF genes, but may be distal",
      x = if (has_log2fc) "Peak accessibility Log2FC" else "Peak-to-gene correlation",
      y = "ZNF gene",
      size = "|Correlation|",
      color = "-log10(link FDR)"
    )
  
  ggsave(
    file.path(out_dir, "Figure_02_ZNF_peak2gene_dotplot.pdf"),
    p_before_dot,
    width = 10,
    height = 7
  )
}

p_before_bar <- znf_gene_summary %>%
  filter(geneName %in% top_znf_genes) %>%
  mutate(geneName = reorder(geneName, sig_rows)) %>%
  ggplot(aes(x = geneName, y = sig_rows)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top ZNF genes by significant peak2gene rows",
    x = "ZNF gene",
    y = "Number of significant rows"
  )

ggsave(
  file.path(out_dir, "Figure_03_ZNF_peak2gene_top_gene_barplot.pdf"),
  p_before_bar,
  width = 8,
  height = 7
)

# =========================================================
# 6. FILTER 3 = KEEP ONLY PHYSICALLY LOCAL PEAKS (ZNF GENE BODY ± 2 KB)
# =========================================================
# This is the actual locus filter.
# Here we ask: among ZNF-linked peaks, which ones overlap the ZNF gene body ± flank_bp?

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

znf_coords <- getBM(
  attributes = c(
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype"
  ),
  filters = "hgnc_symbol",
  values  = unique(znf_p2g_all$geneName),
  mart    = mart
) %>%
  filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
  distinct(hgnc_symbol, chromosome_name, start_position, end_position, .keep_all = TRUE) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name))

cat("ZNF genes retrieved from Ensembl:", nrow(znf_coords), "\n")

write.csv(
  znf_coords,
  file.path(out_dir, "06_ZNF_ensembl_coordinates.csv"),
  row.names = FALSE
)

znf_gr <- GRanges(
  seqnames = znf_coords$chromosome_name,
  ranges = IRanges(
    start = pmax(1, znf_coords$start_position - flank_bp),
    end   = znf_coords$end_position + flank_bp
  ),
  strand = "*",
  geneName = znf_coords$hgnc_symbol
)

# Deduplicate peak coordinates after ranking rows so we retain
# the most significant peak-to-gene link per unique peak.
znf_p2g_unique_peaks <- znf_p2g_all %>%
  rank_peak_links() %>%
  distinct(seqnames, start, end, .keep_all = TRUE)

peaks_gr <- GRanges(
  seqnames = znf_p2g_unique_peaks$seqnames,
  ranges = IRanges(
    start = znf_p2g_unique_peaks$start,
    end   = znf_p2g_unique_peaks$end
  )
)

hits <- findOverlaps(peaks_gr, znf_gr)

znf_local_peaks <- znf_p2g_unique_peaks[queryHits(hits), ] %>%
  mutate(
    locus_gene = znf_gr$geneName[subjectHits(hits)],
    overlaps_local_window = TRUE,
    peak_DA = if (has_log2fc && has_peak_fdr) {
      !is.na(Log2FC) & abs(Log2FC) > log2fc_cutoff & !is.na(FDR.x) & FDR.x < peak_fdr_cutoff
    } else {
      NA
    },
    peak_open = if (has_log2fc && has_peak_fdr) {
      !is.na(Log2FC) & Log2FC > log2fc_cutoff & !is.na(FDR.x) & FDR.x < peak_fdr_cutoff
    } else {
      NA
    },
    peak_closed = if (has_log2fc && has_peak_fdr) {
      !is.na(Log2FC) & Log2FC < -log2fc_cutoff & !is.na(FDR.x) & FDR.x < peak_fdr_cutoff
    } else {
      NA
    }
  )

cat("Unique ZNF-linked peaks overlapping ZNF locus ±2 kb:", nrow(znf_local_peaks), "\n")

write.csv(
  znf_local_peaks,
  file.path(out_dir, "07_ZNF_local_peaks_geneBody_plusminus2kb.csv"),
  row.names = FALSE
)

# =========================================================
# 7. SUMMARIZE THE LOCAL ±2 KB VIEW
# =========================================================
znf_local_counts <- data.frame(
  step = c(
    "Unique ZNF-linked peaks",
    "Unique ZNF-linked peaks overlapping ZNF locus ±2 kb",
    "DA local peaks",
    "Open local peaks",
    "Closed local peaks"
  ),
  count = c(
    nrow(znf_p2g_unique_peaks),
    nrow(znf_local_peaks),
    sum(znf_local_peaks$peak_DA, na.rm = TRUE),
    sum(znf_local_peaks$peak_open, na.rm = TRUE),
    sum(znf_local_peaks$peak_closed, na.rm = TRUE)
  )
)

write.csv(
  znf_local_counts,
  file.path(out_dir, "08_ZNF_local_counts.csv"),
  row.names = FALSE
)

znf_local_summary <- znf_local_peaks %>%
  group_by(locus_gene) %>%
  summarise(
    n_peaks_total = n(),
    n_DA = sum(peak_DA, na.rm = TRUE),
    n_open = sum(peak_open, na.rm = TRUE),
    n_closed = sum(peak_closed, na.rm = TRUE),
    pct_DA = ifelse(n_peaks_total > 0, 100 * n_DA / n_peaks_total, NA),
    mean_Log2FC = if (has_log2fc) mean(Log2FC, na.rm = TRUE) else NA_real_,
    min_peak_FDR = if (has_peak_fdr) min(FDR.x, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(desc(n_DA), desc(mean_Log2FC))

write.csv(
  znf_local_summary,
  file.path(out_dir, "09_ZNF_local_locus_summary.csv"),
  row.names = FALSE
)

top_local_genes <- znf_local_summary %>%
  filter(n_DA > 0) %>%
  slice_head(n = top_n) %>%
  pull(locus_gene)

# =========================================================
# 8. PLOT B = FILTERED LOCAL ±2 KB VIEW
# =========================================================
p_counts_after <- ggplot(znf_local_counts, aes(x = step, y = count, fill = step)) +
  geom_col(show.legend = FALSE) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(
    title = paste0("Local ZNF locus view (gene body ± ", flank_bp / 1000, " kb)"),
    subtitle = "Only peaks physically overlapping the ZNF locus window",
    x = "",
    y = "Number of peaks"
  )

ggsave(
  file.path(out_dir, "Figure_04_ZNF_local_counts.pdf"),
  p_counts_after,
  width = 9,
  height = 5
)

plot_after_df <- znf_local_peaks %>%
  filter(locus_gene %in% top_local_genes)

if (nrow(plot_after_df) > 0 && has_log2fc && has_peak_fdr) {
  plot_after_df <- plot_after_df %>%
    mutate(
      locus_gene = factor(locus_gene, levels = rev(top_local_genes)),
      neglog10_peakFDR = -log10(FDR.x)
    )
  
  p_after_dot <- ggplot(
    plot_after_df,
    aes(
      x = Log2FC,
      y = locus_gene,
      color = neglog10_peakFDR
    )
  ) +
    geom_point(alpha = 0.85, size = 2.5) +
    theme_bw(base_size = 12) +
    labs(
      title = "Local ZNF peaks after ±2 kb filter",
      subtitle = "These peaks physically overlap the ZNF gene body ±2 kb",
      x = "Peak accessibility Log2FC",
      y = "ZNF locus",
      color = "-log10(peak FDR)"
    )
  
  ggsave(
    file.path(out_dir, "Figure_05_ZNF_local_dotplot.pdf"),
    p_after_dot,
    width = 10,
    height = 7
  )
}

p_after_bar <- znf_local_summary %>%
  filter(locus_gene %in% top_local_genes) %>%
  mutate(locus_gene = reorder(locus_gene, n_DA)) %>%
  ggplot(aes(x = locus_gene, y = n_DA)) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    title = "Top local ZNF loci by number of DA peaks",
    x = "ZNF locus",
    y = "Number of DA local peaks"
  )

ggsave(
  file.path(out_dir, "Figure_06_ZNF_local_top_locus_barplot.pdf"),
  p_after_bar,
  width = 8,
  height = 7
)

# =========================================================
# 9. DIRECT BEFORE / AFTER COMPARISON
# =========================================================
comparison_counts <- data.frame(
  stage = c(
    "All ZNF-linked peak2gene rows",
    "Significant ZNF-linked peak2gene rows",
    "Unique ZNF-linked peaks",
    "Local unique peaks (±2 kb)",
    "DA local peaks"
  ),
  count = c(
    nrow(znf_p2g_all),
    nrow(znf_p2g_sig),
    nrow(znf_p2g_unique_peaks),
    nrow(znf_local_peaks),
    sum(znf_local_peaks$peak_DA, na.rm = TRUE)
  )
)

write.csv(
  comparison_counts,
  file.path(out_dir, "10_before_after_comparison_counts.csv"),
  row.names = FALSE
)

p_compare <- ggplot(comparison_counts, aes(x = stage, y = count, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = count), vjust = -0.4, size = 4) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(
    title = "Filter progression from peak2gene ZNF links to local ±2 kb ZNF peaks",
    subtitle = "This makes the difference between linked peaks and physically local peaks explicit",
    x = "",
    y = "Count"
  )

ggsave(
  file.path(out_dir, "Figure_07_filter_progression.pdf"),
  p_compare,
  width = 10,
  height = 5
)

combined_panel <- (p_counts_before | p_counts_after) / p_compare

ggsave(
  file.path(out_dir, "Figure_08_combined_before_after_panel.pdf"),
  combined_panel,
  width = 12,
  height = 10
)

# =========================================================
# 10. LOCAL VS NON-LOCAL AMONG SIGNIFICANT UNIQUE PEAKS
# =========================================================
# This shows the important distinction:
# linked to ZNF gene versus physically local to ZNF locus.
# Duplicate peak coordinates are ranked so the most significant
# peak-to-gene link per unique peak is retained.

sig_unique_peaks <- znf_p2g_sig %>%
  rank_peak_links() %>%
  distinct(seqnames, start, end, .keep_all = TRUE) %>%
  mutate(peak_id = paste(seqnames, start, end, sep = "_"))

znf_local_peaks <- znf_local_peaks %>%
  mutate(peak_id = paste(seqnames, start, end, sep = "_"))

sig_unique_peaks_annot <- sig_unique_peaks %>%
  mutate(
    locus_category = ifelse(
      peak_id %in% znf_local_peaks$peak_id,
      "Local to ZNF locus ±2 kb",
      "Non-local ZNF-linked peak"
    ),
    gene = geneName,
    neglog10_linkFDR = -log10(FDR.y)
  )

write.csv(
  sig_unique_peaks_annot,
  file.path(out_dir, "11_significant_unique_peaks_local_vs_nonlocal.csv"),
  row.names = FALSE
)

local_nonlocal_counts <- sig_unique_peaks_annot %>%
  count(locus_category)

write.csv(
  local_nonlocal_counts,
  file.path(out_dir, "12_local_vs_nonlocal_counts.csv"),
  row.names = FALSE
)

p_local_nonlocal_bar <- ggplot(local_nonlocal_counts, aes(x = locus_category, y = n, fill = locus_category)) +
  geom_col(show.legend = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "Among significant ZNF-linked peaks, which are physically local?",
    subtitle = "Local = overlaps ZNF gene body ±2 kb",
    x = "",
    y = "Number of unique significant peaks"
  )

ggsave(
  file.path(out_dir, "Figure_09_local_vs_nonlocal_barplot.pdf"),
  p_local_nonlocal_bar,
  width = 7,
  height = 5
)

gene_order2 <- sig_unique_peaks_annot %>%
  group_by(gene) %>%
  summarise(mean_Log2FC = mean(Log2FC, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_Log2FC) %>%
  pull(gene)

sig_unique_peaks_annot <- sig_unique_peaks_annot %>%
  mutate(gene = factor(gene, levels = gene_order2))

p_local_nonlocal_dot <- ggplot(
  sig_unique_peaks_annot,
  aes(x = Log2FC, y = gene)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(
    aes(size = abs(Log2FC), color = neglog10_linkFDR, shape = locus_category),
    alpha = 0.9
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Significant ZNF-linked peaks: local vs non-local",
    subtitle = "This separates distal peak2gene links from physically local ZNF peaks",
    x = "Peak accessibility Log2FC",
    y = "ZNF gene",
    size = "|Log2FC|",
    color = "-log10(link FDR)",
    shape = ""
  )

ggsave(
  file.path(out_dir, "Figure_10_local_vs_nonlocal_dotplot.pdf"),
  p_local_nonlocal_dot,
  width = 11,
  height = 12
)

# =========================================================
# 11. EXPORT A SIMPLE SUMMARY TEXT FILE
# =========================================================
summary_lines <- c(
  paste("Input peak2gene rows:", nrow(df)),
  paste("ZNF-linked peak2gene rows:", nrow(znf_p2g_all)),
  paste("Significant ZNF-linked peak2gene rows:", nrow(znf_p2g_sig)),
  paste("Unique ZNF-linked peaks:", nrow(znf_p2g_unique_peaks)),
  paste("Unique peaks local to ZNF locus ±2 kb:", nrow(znf_local_peaks)),
  paste("DA local peaks:", sum(znf_local_peaks$peak_DA, na.rm = TRUE)),
  paste("Unique ZNF genes in peak2gene links:", length(unique(znf_p2g_all$geneName))),
  paste("Unique ZNF loci with local overlapping peaks:", length(unique(znf_local_peaks$locus_gene)))
)

writeLines(
  summary_lines,
  con = file.path(out_dir, "13_summary.txt")
)

cat("\nDone.\n")
cat("Outputs written to:", out_dir, "\n")

# =========================================================
# 12. EXPORT LOCAL SIGNIFICANT PEAKS AS BED / GFF
# =========================================================
local_sig_peaks <- sig_unique_peaks_annot %>%
  filter(locus_category == "Local to ZNF locus ±2 kb")

# Create BED format
bed_df <- local_sig_peaks %>%
  mutate(
    chrom = seqnames,
    chromStart = start - 1,   # BED is 0-based
    chromEnd = end,
    name = paste0(gene, "_", round(Log2FC, 2)),
    score = 1000,
    strand = "."
  ) %>%
  dplyr::select(chrom, chromStart, chromEnd, name, score, strand)

write.table(
  bed_df,
  file = file.path(out_dir, "ZNF_local_21_peaks.bed"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

gff_df <- local_sig_peaks %>%
  mutate(
    seqid = seqnames,
    source = "ATAC",
    type = "peak",
    start = start,
    end = end,
    score = ".",
    strand = ".",
    phase = ".",
    attributes = paste0("ID=", gene, ";Log2FC=", round(Log2FC, 2))
  ) %>%
  dplyr::select(seqid, source, type, start, end, score, strand, phase, attributes)

write.table(
  gff_df,
  file = file.path(out_dir, "ZNF_local_21_peaks.gff"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



# =========================================================
# EXPORT TOP HIT PER PEAK (BEST LINK PER REGION)
# =========================================================

top_hits_per_peak <- sig_unique_peaks %>%
  arrange(FDR.y) %>%   # ensure best links are prioritized
  mutate(
    peak_id = paste(seqnames, start, end, sep = "_")
  )

write.csv(
  top_hits_per_peak,
  file.path(out_dir, "TOP_HITS_per_peak_best_link.csv"),
  row.names = FALSE
)
