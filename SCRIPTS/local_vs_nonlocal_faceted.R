library(tidyverse)
library(biomaRt)

files <- list.files("PEAKv3", pattern = "\\.csv$", full.names = TRUE)

all_znf_genes <- files %>%
  map_dfr(~ read.csv(.x) %>% filter(grepl("^ZNF", geneName))) %>%
  pull(geneName) %>%
  unique()

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

znf_coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = all_znf_genes,
  mart = mart
) %>%
  filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
  distinct()

write.csv(znf_coords, "ZNF_ensembl_coordinates_cached.csv", row.names = FALSE)


# library(biomaRt)
#
# mart <- useMart(
#   biomart = "ENSEMBL_MART_ENSEMBL",
#   dataset = "hsapiens_gene_ensembl",
#   host = "https://www.ensembl.org"
# )
#
# znf_coords <- getBM(
#   attributes = c(
#     "hgnc_symbol",
#     "chromosome_name",
#     "start_position",
#     "end_position",
#     "strand",
#     "gene_biotype"
#   ),
#   filters = "hgnc_symbol",
#   values = all_znf_genes,
#   mart = mart
# )

#######################################


process_file <- function(file, znf_coords) {

  df <- read.csv(file, stringsAsFactors = FALSE)

  condition <- basename(file) %>%
    stringr::str_remove("v3_significant_peaks_") %>%
    stringr::str_remove("\\.csv")

  df <- df %>%
    mutate(seqnames = as.character(seqnames))

  znf_coords <- znf_coords %>%
    mutate(
      chromosome_name = ifelse(
        grepl("^chr", chromosome_name),
        chromosome_name,
        paste0("chr", chromosome_name)
      )
    )

  if (!any(grepl("^chr", df$seqnames)) && any(grepl("^chr", znf_coords$chromosome_name))) {
    df <- df %>%
      mutate(seqnames = paste0("chr", seqnames))
  }

  znf <- df %>%
    filter(grepl("^ZNF", geneName)) %>%
    mutate(
      is_sig =
        !is.na(Correlation) &
        abs(Correlation) > 0.2 &
        !is.na(FDR.y) &
        FDR.y < 0.05 &
        !is.na(Log2FC) &
        abs(Log2FC) > 0.5 &
        !is.na(FDR.x) &
        FDR.x < 0.05
    )

  znf_sig_unique <- znf %>%
    filter(is_sig) %>%
    arrange(FDR.y, desc(abs(Correlation)), desc(abs(Log2FC))) %>%
    distinct(seqnames, start, end, .keep_all = TRUE) %>%
    mutate(peak_id = paste(seqnames, start, end, sep = "_"))

  if (nrow(znf_sig_unique) == 0) {
    return(tibble(
      seqnames = NA_character_,
      start = NA_real_,
      end = NA_real_,
      geneName = "No significant peaks",
      Log2FC = 0,
      Correlation = NA_real_,
      FDR.x = NA_real_,
      FDR.y = NA_real_,
      peak_id = NA_character_,
      locus_category = "No significant ZNF-linked peaks",
      neglog10_linkFDR = NA_real_,
      condition = condition,
      is_placeholder = TRUE
    ))
  }

  znf_gr <- GenomicRanges::GRanges(
    seqnames = znf_coords$chromosome_name,
    ranges = IRanges::IRanges(
      start = pmax(1, znf_coords$start_position - 2000),
      end = znf_coords$end_position + 2000
    )
  )

  peaks_gr <- GenomicRanges::GRanges(
    seqnames = znf_sig_unique$seqnames,
    ranges = IRanges::IRanges(
      start = znf_sig_unique$start,
      end = znf_sig_unique$end
    )
  )

  hits <- GenomicRanges::findOverlaps(peaks_gr, znf_gr)

  local_peak_ids <- znf_sig_unique$peak_id[queryHits(hits)]

  znf_sig_unique %>%
    mutate(
      locus_category = ifelse(
        peak_id %in% local_peak_ids,
        "Local to ZNF locus ±2 kb",
        "Non-local ZNF-linked peak"
      ),
      locus_category = as.character(locus_category),
      neglog10_linkFDR = -log10(FDR.y),
      condition = condition,
      is_placeholder = FALSE
    )
}








##############################
# =========================================================
# FINAL PLOT SECTION
# =========================================================

all_data <- map_dfr(files, process_file, znf_coords = znf_coords)

write.csv(
  all_data,
  "ALL_local_vs_nonlocal_results_with_empty_panels.csv",
  row.names = FALSE
)

gene_order <- all_data %>%
  filter(is_placeholder == FALSE) %>%
  group_by(geneName) %>%
  summarise(mean_Log2FC = mean(Log2FC, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_Log2FC) %>%
  pull(geneName)

all_data <- all_data %>%
  mutate(
    geneName = ifelse(is_placeholder, "No significant peaks", geneName),
    geneName = factor(geneName, levels = c(gene_order, "No significant peaks")),
    condition = factor(
      condition,
      levels = c(
        "BG_CD14_Mono",
        "BG_CD4_TCM",
        "MDP_CD14_Mono",
        "MDP_CD4_TCM",
        "SYKi_CD14_Mono",
        "SYKi_CD4_TCM"
      )
    )
  )

p_final <- ggplot() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(
    data = all_data %>% filter(is_placeholder == FALSE),
    aes(
      x = Log2FC,
      y = geneName,
      size = abs(Log2FC),
      color = neglog10_linkFDR,
      shape = locus_category
    ),
    alpha = 0.85
  ) +
  geom_text(
    data = all_data %>% filter(is_placeholder == TRUE),
    aes(
      x = 0,
      y = geneName,
      label = "No significant peaks"
    ),
    size = 4
  ) +
  facet_wrap(~ condition, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Significant ZNF-linked peaks: local vs non-local",
    subtitle = "Empty panels indicate no ZNF-linked peaks passed correlation, FDR, and Log2FC thresholds",
    x = "Peak accessibility Log2FC",
    y = "ZNF gene",
    size = "|Log2FC|",
    color = "-log10(link FDR)",
    shape = ""
  )

ggsave(
  "FINAL_local_vs_nonlocal_faceted_with_empty_panels.pdf",
  p_final,
  width = 16,
  height = 12,
  limitsize = FALSE
)

