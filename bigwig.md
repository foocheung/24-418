# BIGWIG Generation Pipeline (ATAC-seq and scRNA-seq)

This repository documents the generation of BIGWIG tracks from ATAC-seq (ArchR) and scRNA-seq (10x Genomics) data for downstream visualization and comparative analysis.

---

## Overview

Two types of BIGWIG tracks were generated:

1. **ATAC-seq BIGWIGs**
   - Generated using ArchR
   - Represent chromatin accessibility
   - Grouped by cell type and treatment

2. **scRNA-seq BIGWIGs**
   - Generated from 10x BAM files
   - Represent expression coverage
   - Split by cell type using Seurat annotations

Scripts used can be found here: https://github.com/foocheung/24-418/tree/main/SCRIPTS/BigWigs

---

## ATAC-seq BIGWIG Generation

### Input
- ArchRProject objects (`.rds`)
- Peak matrix and fragment files generated within ArchR

### Processing Steps

1. Load ArchR project:
```r
proj <- readRDS("proj_scrnaseq_unconstrained.rds")
```

2. Annotate treatment:

```r
proj$Treatment <- ifelse(grepl("^PMA_POS", proj$Sample), "PMA_POS", "PMA_NEG")
```

3. Define grouping:

```r
proj$group_bw <- paste(proj$Clusters2, proj$Treatment, sep = "_")
```

4. Generate BIGWIG tracks:

```r
getGroupBW(
  ArchRProj = proj,
  groupBy = "group_bw",
  normMethod = "ReadsInTSS",
  tileSize = 25
)
```

### Output

* Location:

```
ArchRProject/GroupBigWigs/group_bw/
```

* File format:

```
<celltype>_<PMA_STATUS>-TileSize-25-normMethod-ReadsInTSS-ArchR.bw
```

### Notes

* Normalization: **ReadsInTSS**
* Resolution: **25 bp tiles**
* Represents **pseudo-bulk chromatin accessibility**

---

## scRNA-seq BIGWIG Generation

### Input

* 10x Genomics BAM files (`possorted_genome_bam.bam`)
* Seurat object with cell annotations:

```
data2_ref.rds
```

### Key Metadata

* `predicted.celltype.l2` → cell type labels
* `Lane` → maps cells to GEX1–GEX8 BAMs

---

### Step 1: Extract cell barcodes

Seurat cell names:

```
AAACCCAAGACGGTTG-2_24_418
```

Converted to BAM-compatible barcodes:

```r
barcode_merged <- sub("_.*$", "", cell)
barcode <- sub("-[0-9]+$", "-1", barcode_merged)
```

This ensures matching with BAM tags:

```
CB:Z:AAACCCAAGACGGTTG-1
```

---

### Step 2: Split BAM by cell type

For each GEX lane:

* Subset reads using barcode lists
* Match against `CB:Z:` tag

---

### Step 3: Generate BIGWIGs

Using `bamCoverage`:

```bash
bamCoverage \
  --bam subset.bam \
  --outFileName output.bw \
  --outFileFormat bigwig \
  --binSize 25 \
  --normalizeUsing CPM \
  --numberOfProcessors 8 \
  --minMappingQuality 30 \
  --skipNonCoveredRegions
```

### Output

* One BIGWIG per:

  * GEX lane
  * Cell type

Example:

```
GEX1/CD4_TCM.bw
GEX2/CD14_Mono.bw
```

---

## Normalization

| Data Type | Method     | Purpose                                         |
| --------- | ---------- | ----------------------------------------------- |
| ATAC-seq  | ReadsInTSS | Controls for TSS enrichment and library quality |
| scRNA-seq | CPM        | Normalizes for sequencing depth                 |

---

## Biological Notes

* BG-trained immunity effects were strongest in **monocytes**
* Significant chromatin accessibility changes were primarily observed in:

  * **CD14 Monocytes**
* Other cell types showed weaker signals and often did not pass significance thresholds
* This motivated downstream **motif-based analyses**

---

## Directory Structure

```
BIGWIGS/
├── ATAC/
│   ├── BG/
│   ├── SYKI/
│   ├── UNTRAINED/
│   ├── MDP/
├── SCRNA/
│   ├── GEX1/
│   ├── GEX2/
│   ├── GEX3/
│   ├── GEX4/
│   ├── GEX5/
│   ├── GEX6/
│   ├── GEX7/
│   ├── GEX8/
```

---

## Reproducibility

### Requirements

* R (Seurat, ArchR, tidyverse)
* samtools
* deepTools (`bamCoverage`)

### Key Steps

1. Generate ATAC BIGWIGs via ArchR
2. Extract barcodes from Seurat
3. Split BAM by cell type
4. Run `bamCoverage` per subset
5. Organize outputs

---

## Logs and Provenance

Full file paths and logs for BIGWIG generation are documented here:


---

## Summary

* ATAC BIGWIGs → chromatin accessibility (ArchR, ReadsInTSS)
* scRNA BIGWIGs → expression coverage (bamCoverage, CPM)
* Cell-type–specific tracks generated across all datasets
* Monocytes show strongest trained immunity signal

---

