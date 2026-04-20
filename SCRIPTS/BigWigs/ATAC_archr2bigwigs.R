library(ArchR)
library(tidyverse)

addArchRThreads(threads = 8)

# ------------------------------------------------------------------------------
# Helper function to add treatment/group labels and generate group bigWigs
# ------------------------------------------------------------------------------

make_group_bigwigs <- function(
  project_rds,
  threads = 8,
  tile_size = 25,
  norm_method = "ReadsInTSS"
) {
  proj <- readRDS(project_rds)

  proj$Treatment <- ifelse(
    grepl("^PMA_POS", proj$Sample),
    "PMA_POS",
    "PMA_NEG"
  )

  proj$group_bw <- paste(proj$Clusters2, proj$Treatment, sep = "_")

  addArchRThreads(threads = threads)

  getGroupBW(
    ArchRProj = proj,
    groupBy = "group_bw",
    normMethod = norm_method,
    tileSize = tile_size
  )

  invisible(proj)
}

# ------------------------------------------------------------------------------
# Project paths
# ------------------------------------------------------------------------------

project_paths <- list(
  bg   = "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/a7/961ca94603d713bbcbc755966180e8/proj_scrnaseq_unconstrained.rds",
  ski  = "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/d2/3855bf25254f5e25ece914412c904c/proj_scrnaseq_unconstrained.rds",
  syki = "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/d2/3855bf25254f5e25ece914412c904c/proj_scrnaseq_unconstrained.rds",
  un   = "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/d2/3855bf25254f5e25ece914412c904c/proj_scrnaseq_unconstrained.rds",
  mdp  = "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/a7/05f2f92f55e67bc78318f4c3065114/proj_scrnaseq_unconstrained.rds"
)

# ------------------------------------------------------------------------------
# Run group bigWig generation for each project
# ------------------------------------------------------------------------------

project_list <- purrr::imap(
  project_paths,
  ~ {
    message("Processing: ", .y)
    make_group_bigwigs(project_rds = .x)
  }
)

# ------------------------------------------------------------------------------
# Example output directory
# ------------------------------------------------------------------------------

# Group bigWigs are typically written under:
# <ArchRProject>/GroupBigWigs/group_bw/
#
# Example:
# /data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/16/49788349988a8ae13aa5568257c52b/ArchRProject/GroupBigWigs/group_bw/
