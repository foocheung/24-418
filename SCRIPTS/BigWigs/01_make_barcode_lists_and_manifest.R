suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

rds_file <- "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/assets/data2_ref.rds"
runroot  <- "/data/chi/PROJECTS/24-418_Fraser_HIV_Reactivation/Bioinformatics/Runs/22MN7CLT3/RUN2"
workroot <- "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe"

a <- readRDS(rds_file)

stopifnot("Lane" %in% colnames(a@meta.data))
stopifnot("predicted.celltype.l2" %in% colnames(a@meta.data))

meta <- a@meta.data
meta$cell <- rownames(meta)

# Example merged cell name:
# AAACCCAAGACGGTTG-2_24_418
# First remove the trailing sample suffix, then convert merged suffix back to -1
meta$barcode_merged <- sub("_.*$", "", meta$cell)
meta$barcode <- sub("-[0-9]+$", "-1", meta$barcode_merged)

meta <- meta %>%
  mutate(
    Lane = as.character(Lane),
    predicted.celltype.l2 = as.character(predicted.celltype.l2)
  ) %>%
  filter(Lane %in% as.character(1:8))

meta$celltype_safe <- meta$predicted.celltype.l2 %>%
  str_replace_all("[[:space:]/]+", "_") %>%
  str_replace_all("[^A-Za-z0-9_\\-]", "")

barcode_root <- file.path(workroot, "barcode_lists_by_GEX")
dir.create(barcode_root, recursive = TRUE, showWarnings = FALSE)

bam_map <- c(
  GEX1 = file.path(runroot, "3config_24_418_GEX1/outs/possorted_genome_bam.bam"),
  GEX2 = file.path(runroot, "3config_24_418_GEX2/outs/possorted_genome_bam.bam"),
  GEX3 = file.path(runroot, "3config_24_418_GEX3/outs/possorted_genome_bam.bam"),
  GEX4 = file.path(runroot, "3config_24_418_GEX4/outs/possorted_genome_bam.bam"),
  GEX5 = file.path(runroot, "3config_24_418_GEX5/outs/possorted_genome_bam.bam"),
  GEX6 = file.path(runroot, "3config_24_418_GEX6/outs/possorted_genome_bam.bam"),
  GEX7 = file.path(runroot, "3config_24_418_GEX7/outs/possorted_genome_bam.bam"),
  GEX8 = file.path(runroot, "3config_24_418_GEX8/outs/possorted_genome_bam.bam")
)

manifest <- vector("list", length = 0)

for (lane_i in sort(unique(meta$Lane))) {
  gex <- paste0("GEX", lane_i)
  gex_dir <- file.path(barcode_root, gex)
  dir.create(gex_dir, recursive = TRUE, showWarnings = FALSE)

  meta_lane <- meta %>% filter(Lane == lane_i)

  split_list <- split(meta_lane$barcode, meta_lane$celltype_safe)

  write.csv(
    meta_lane[, c("cell", "barcode_merged", "barcode", "Lane", "predicted.celltype.l2", "celltype_safe")],
    file = file.path(gex_dir, "cell_barcode_mapping.csv"),
    row.names = FALSE
  )

  for (ct in names(split_list)) {
    bc <- unique(split_list[[ct]])
    bc <- bc[!is.na(bc) & bc != ""]
    if (length(bc) == 0) next

    bc_file <- file.path(gex_dir, paste0(ct, ".txt"))
    writeLines(bc, bc_file)

    outdir <- file.path(runroot, gex)
    tmpdir <- file.path(outdir, "tmp_bams")
    logdir <- file.path(outdir, "logs")

    manifest[[length(manifest) + 1]] <- data.frame(
      gex = gex,
      lane = lane_i,
      celltype = ct,
      barcode_file = bc_file,
      bam = unname(bam_map[gex]),
      outdir = outdir,
      tmpdir = tmpdir,
      logdir = logdir,
      bw = file.path(outdir, paste0(ct, ".bw")),
      sorted_bam = file.path(tmpdir, paste0(ct, ".sorted.bam")),
      stringsAsFactors = FALSE
    )
  }
}

manifest_df <- dplyr::bind_rows(manifest)

manifest_file <- file.path(workroot, "bigwig_manifest.tsv")
write.table(
  manifest_df,
  file = manifest_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Wrote manifest:", manifest_file, "\n")
cat("Number of jobs:", nrow(manifest_df), "\n")
