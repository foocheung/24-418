## Analysis Scripts Overview

This repository contains scripts for integrated **scRNA-seq and scATAC-seq analysis**, including pseudobulk differential expression, pathway enrichment, transcription factor analysis, and genome-level visualization.

---

### scRNA-seq Analysis

#### `pseudobulk_deltaGSEA_pipeline.R`

Performs pseudobulk differential expression analysis across training conditions (Untrained, SYKi, MDP, BG) and PMA treatment.
Implements a **delta-of-delta log2FC framework** to isolate training-specific responses relative to the untrained baseline, followed by GSEA using BTM pathways.

#### `scRNAseq_GSEA_bubble_plot.R`

Generates bubble plot visualizations of GSEA results.

* Color: Normalized Enrichment Score (NES)
* Size: −log10(adjusted p-value)
* Faceted by cell type
  Produces plots for individual and combined conditions.

#### `scRNAseq_module_gene_feature_violin_plots.R`

Visualizes gene expression and module scores using:

* FeaturePlots (embedding-based)
* Violin plots (distribution summaries)

Supports:

* Module scores (`ModuleScore*`)
* Gene expression (RNA assay)

Includes:

* Quantile-based scaling
* Median annotations
* Optional non-zero filtering for gene expression

---

### ATAC-seq Analysis

#### `ATAC_peak2gene_pathway_enrichment.R`

Performs pathway enrichment using genes linked to ATAC-seq peaks.
Filters peak-gene links based on:

* Accessibility change (Log2FC)
* Correlation
* Link FDR

Runs Reactome and GO enrichment to identify pathways associated with chromatin accessibility.

---

#### `ATAC_TF_motif_family_analysis.R`

Analyzes transcription factor (TF) motif activity and family composition.

* Uses chromVAR-derived TFs
* Maps TFs to families (CisBP database)
* Aggregates TFs **irrespective of direction (up/down)** to assess overall family representation
* Visualizes overlaps using UpSet plots and bar charts

---

#### `ATAC_ZNF_peak2gene_local_analysis.R`

Distinguishes between:

* **Peak2gene-linked peaks** (regulatory associations, potentially distal)
* **Physically local peaks** (overlapping ZNF gene loci ±2 kb)

Includes:

* Ranking of peak-gene links
* Ensembl-based gene coordinate retrieval
* Local vs non-local classification
* Summary statistics and comparative plots

---

#### `ATAC_archr2bigwigs.R`

Exports ArchR ATAC-seq signal tracks to BigWig format for genome browser visualization.
Used for downstream visualization in IGV or UCSC Genome Browser.

---

### Genome-Level Visualization

#### `ATAC_genome_circos_visualization.R`

Generates genome-wide circos plots showing:

* Binned ATAC accessibility across chromosomes
* Differential peaks
* Optional HIV integration overlays

Used for global comparison of chromatin accessibility patterns.

---

#### `HIV_ATAC_genome_browser_tracks.R`

Creates genome browser-style plots of:

* HIV gene annotations
* ATAC-seq accessibility tracks across conditions

Aligns all tracks to the HIV genome for direct visualization of accessibility patterns.

---

#### `Peak–Gene_linking_and_HIV_Integration_Overlap_Pipeline.R`

Integrates ATAC peak-to-gene links with HIV integration data.
Identifies overlaps between:

* Accessible chromatin regions
* HIV integration sites
* Gene regulatory regions

---

## Key Concepts Across Scripts

* **Pseudobulk analysis** reduces single-cell noise while preserving biological structure
* **Delta-of-delta framework** isolates training-specific transcriptional responses
* **Peak2gene linking** captures regulatory relationships (not necessarily physical proximity)
* **Local locus filtering** distinguishes true chromatin accessibility at gene loci
* **TF family aggregation** focuses on regulatory programs rather than individual TF direction
* **Visualization scripts** are descriptive and do not perform statistical testing

---

