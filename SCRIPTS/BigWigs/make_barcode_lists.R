suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

rds_file <- "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/assets/data2_ref.rds"
a <- readRDS(rds_file)

stopifnot("Lane" %in% colnames(a@meta.data))
stopifnot("predicted.celltype.l2" %in% colnames(a@meta.data))

meta <- a@meta.data
meta$cell <- rownames(meta)

# For names like: AAACCCAAGACGGTTG-1_24_418
meta$barcode <- sub("_.*$", "", meta$cell)

meta <- meta %>%
  mutate(
    Lane = as.character(Lane),
    predicted.celltype.l2 = as.character(predicted.celltype.l2)
  ) %>%
  filter(Lane %in% as.character(1:8))

meta$celltype_safe <- meta$predicted.celltype.l2 %>%
  str_replace_all("[[:space:]/]+", "_") %>%
  str_replace_all("[^A-Za-z0-9_\\-]", "")

out_root <- "barcode_lists_by_GEX"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

for (lane_i in sort(unique(meta$Lane))) {
  gex_dir <- file.path(out_root, paste0("GEX", lane_i))
  dir.create(gex_dir, showWarnings = FALSE, recursive = TRUE)

  meta_lane <- meta %>% filter(Lane == lane_i)

  split_list <- split(meta_lane$barcode, meta_lane$celltype_safe)

  for (ct in names(split_list)) {
    bc <- unique(split_list[[ct]])
    bc <- bc[!is.na(bc) & bc != ""]
    writeLines(bc, con = file.path(gex_dir, paste0(ct, ".txt")))
  }

  write.csv(
    meta_lane[, c("cell", "barcode", "Lane", "predicted.celltype.l2", "celltype_safe")],
    file = file.path(gex_dir, "cell_barcode_mapping.csv"),
    row.names = FALSE
  )
}

cat("Done writing barcode lists.\n")

