## Methods

### Processing scRNA-seq

Raw 10x Genomics scRNA-seq data were imported into Seurat objects per sequencing lane, with sample identities assigned using Demuxalot. Putative doublets were detected using DoubletFinder and removed prior to downstream processing. Lane-level objects were normalized, variable features identified, and merged into a single Seurat object for integration. Batch correction and integration were performed with Seurat’s CCA, RPCA and Harmony. Cell type annotation was performed using reference-based mapping with the Monaco immune reference dataset, followed by manual curation based on canonical marker gene expression to refine major immune populations and memory subsets. The final annotated object was used for all downstream analyses.

---

### Processing of scATAC-seq data

Single-cell ATAC-seq data were processed using the ArchR framework. Raw fragment files were imported into ArchR projects, and quality control filtering was applied to remove low-quality nuclei based on transcription start site (TSS) enrichment and fragment counts. Putative doublets were identified and removed using ArchR’s doublet detection algorithm.

Dimensionality reduction was performed using iterative latent semantic indexing (LSI), followed by graph-based clustering and UMAP visualization. Peaks were called using MACS2 within ArchR on pseudobulk profiles generated from clustered cells.

To assign cell identities, gene activity scores were computed from chromatin accessibility profiles using ArchR. Cell type labels were transferred from the matched scRNA-seq dataset to scATAC-seq cells using ArchR’s integration framework, enabling consistent annotation across modalities. Transferred labels were validated using accessibility at canonical marker genes.

Peak-to-gene linkages were computed using correlation-based methods implemented in ArchR to associate distal regulatory elements with putative target genes. Motif enrichment analysis was performed using ArchR’s motif annotation framework, and transcription factor activity was inferred using chromVAR deviation scores.

---

### Pathway Enrichment Analysis (scRNA-seq)

Single-cell RNA-seq data were aggregated into pseudobulk expression profiles using Seurat’s AggregateExpression() function. Aggregation was performed by treatment status, cell type (monaco_ann1 annotation), subject (sCALL), and sequencing lane to preserve biological and technical structure in the summarized profiles.

Differential expression analysis was performed on pseudobulk profiles comparing PMA-treated versus untreated samples within each condition. For each trained condition (SYKi, MDP, BG), treatment-induced changes were compared with those observed in the untrained baseline by computing a delta-of-delta log2 fold change:

Δlog2FC = log2FC(PMA-treated vs. untreated in trained condition) − log2FC(PMA-treated vs. untreated in untrained condition).

This approach identifies genes whose stimulation response is increased or decreased in trained conditions relative to the untrained baseline. Genes detected in both the trained and the untrained comparisons were retained for ranking without additional significance filtering.

Ranked gene lists were analyzed using Gene Set Enrichment Analysis (GSEA) implemented in the clusterProfiler package against Blood Transcriptional Modules (BTM). P-values were adjusted using the Benjamini–Hochberg method.

Pathway enrichment results were generated for each cell type and training condition and summarized to compare condition-specific enrichment patterns across SYKi, MDP, BG, and the untrained baseline.

Filtered GSEA results were visualized as bubble plots using ggplot2. Pathways with adjusted p-values less than 0.05 were retained, and entries containing “TBA” in the pathway ID were excluded. Bubble plots were generated separately for MDP/BG, SYKi, and all conditions combined. In each plot, the x-axis represents training condition, the y-axis represents enriched pathway, bubble color represents normalized enrichment score (NES), and bubble size represents statistical significance as -log10(p.adjust). Results were faceted by cell type using the monaco_ann1 annotation, with cell types displayed in a manually defined order. Plots were exported as both PNG and SVG files for manuscript figure generation.

---

### ATAC-seq Associated Gene Pathway Analysis

Differential chromatin accessibility analysis revealed that statistically significant peaks passing multiple testing correction (FDR-adjusted threshold) were most robustly detected in β-glucan (BG)-trained CD14⁺ monocytes, whereas other training conditions and cell types showed fewer peaks meeting stringent significance thresholds. To enable comparative pathway-level interpretation across all training conditions, a peak-to-gene linkage–based enrichment approach was applied using correlated accessibility changes, allowing inclusion of biologically relevant regulatory signals that may not reach peak-level significance in all conditions. This approach is particularly suited to single-cell ATAC-seq data, where peak-level statistical power can be limited due to data sparsity.

Pathway enrichment analysis was performed on genes linked to differentially accessible ATAC-seq peaks using the clusterProfiler and ReactomePA packages. Peak-to-gene associations were derived from ArchR-generated correlation-based linkages between chromatin accessibility and gene activity scores. Peak–gene pairs were filtered using thresholds on differential accessibility (log2 fold change) and peak–gene correlation strength to retain high-confidence regulatory associations. Gene symbols were converted to ENTREZ identifiers using the org.Hs.eg.db annotation package.

Reactome and Gene Ontology (GO) enrichment analyses were performed on the resulting gene sets. Statistical significance was assessed using Benjamini–Hochberg correction, and pathways with adjusted p-values below the significance threshold were retained. Enrichment results were visualized using bar plots (Reactome) and dot plots (GO). Analyses were conducted separately for each training condition and cell type.

---

### ATAC-seq Peak-to-Gene and Locus-Level Accessibility Analysis of ZNF Transcription Factors

Peak-to-gene linkage analysis was performed using ArchR in β-glucan (BG)-trained versus untrained conditions. Peak–gene links were filtered based on correlation (>0.2), linkage FDR (<0.05), peak-level FDR (<0.05), and |Log2FC| (>0.5).

ZNF-associated links were identified by selecting genes with the “ZNF” prefix. Peak-to-gene links were ranked and the most confident gene assignment per peak was retained.

ZNF gene coordinates were retrieved using biomaRt, extended by ±2 kb, and intersected with ATAC peaks using GenomicRanges. Peaks were classified as locus-proximal or distal.

---

### Transcription Factor Motif Activity and Family Analysis

TF motif activity was quantified using chromVAR deviation scores implemented in ArchR. Differences in motif activity were assessed using Wilcoxon rank-sum tests with FDR correction.

TFs were mapped to families using the CisBP database. TFs were aggregated irrespective of direction of change to assess overall family representation.

Overlap across conditions was visualized using UpSet plots, with family composition summarized using bar and stacked plots.

---

### Feature Visualization of Gene Expression and Module Scores

Gene expression and module scores were visualized using FeaturePlot and violin plots. Cells were grouped by training condition and treatment status.

Expression values were quantile-scaled. Violin plots included median annotations and optional non-zero filtering.

---

### Genome-wide Circos Visualization of ATAC Accessibility and HIV Integration

Genome-wide ATAC accessibility and HIV integration were visualized using circos plots. Accessibility was shown as radial bars, integration sites as points, and peaks as line tracks.

---

### Genome Browser Visualization

HIV genome browser-style plots were generated by mapping ATAC peaks and HIV gene annotations. Tracks were aligned across conditions and exported as PDF/SVG.

---

### Chromatin Peak–HIV Integration Overlap Framework

Peak–gene links and HIV integration sites were intersected to identify integration within accessible chromatin. Peak–gene associations were stratified by correlation and analyzed using pathway enrichment.

---

### ZNF Peak-to-Gene and Local Locus Accessibility Analysis

Peak-to-gene links were filtered for ZNF genes. Peaks were classified as local (±2 kb overlap) or distal. Summary statistics and visualization were generated.
