#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript run_group_bw.R <project_rds> <threads>")
}

project_rds <- args[1]
threads <- as.integer(args[2])

message("Reading project: ", project_rds)
proj <- readRDS(project_rds)

# Basic checks
cell_cols <- colnames(getCellColData(proj))
message("Available cellColData columns: ", paste(cell_cols, collapse = ", "))

if (!"Sample" %in% cell_cols) stop("Sample column not found in ArchR project.")
if (!"Clusters2" %in% cell_cols) stop("Clusters2 column not found in ArchR project.")

# Build grouping
proj$Treatment <- ifelse(grepl("^PMA_POS", proj$Sample), "PMA_POS", "PMA_NEG")
proj$group_bw <- paste(gsub(" ", "_", proj$Clusters2), proj$Treatment, sep = "_")

message("Group counts:")
print(table(proj$group_bw))

addArchRThreads(threads = threads)

# Generate BigWigs
getGroupBW(
  ArchRProj = proj,
  groupBy = "group_bw",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  force = TRUE
)

message("Done: ", project_rds)
message("Output dir: ", outputDirectory(proj))
