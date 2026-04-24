## Methods

### Processing scRNA-seq

Raw 10x Genomics scRNA-seq data were imported into Seurat objects per sequencing lane, with sample identities assigned using Demuxalot. Putative doublets were detected using DoubletFinder and removed prior to downstream processing. Lane-level objects were normalized, variable features identified, and merged into a single Seurat object for integration. Batch correction and integration were performed with Seurat’s CCA, RPCA and Harmony. Cell type annotation was performed using reference-based mapping with the Monaco immune reference dataset, followed by manual curation based on canonical marker gene expression to refine major immune populations and memory subsets. The final annotated object was used for all downstream analyses.

### Processing of scATAC-seq data

Single-cell ATAC-seq data were processed using the ArchR framework. Raw fragment files were imported into ArchR projects, and quality control filtering was applied to remove low-quality nuclei based on transcription start site (TSS) enrichment and fragment counts. Putative doublets were identified and removed using ArchR’s doublet detection algorithm.

Dimensionality reduction was performed using iterative latent semantic indexing (LSI), followed by graph-based clustering and UMAP visualization. Peaks were called using MACS2 within ArchR on pseudobulk profiles generated from clustered cells.

To assign cell identities, gene activity scores were computed from chromatin accessibility profiles using ArchR. Cell type labels were transferred from the matched scRNA-seq dataset to scATAC-seq cells using ArchR’s integration framework, enabling consistent annotation across modalities. Transferred labels were validated using accessibility at canonical marker genes.

Peak-to-gene linkages were computed using correlation-based methods implemented in ArchR to associate distal regulatory elements with putative target genes. Motif enrichment analysis was performed using ArchR’s motif annotation framework, and transcription factor activity was inferred using chromVAR deviation scores.

### Pathway Enrichment Analysis (scRNA-seq)

Single-cell RNA-seq data were aggregated into pseudobulk expression profiles using Seurat’s `AggregateExpression()` function. Aggregation was performed by treatment status, cell type (`monaco_ann1` annotation), subject (`sCALL`), and sequencing lane to preserve biological and technical structure in the summarized profiles.

Differential expression analysis was performed on pseudobulk profiles comparing PMA-treated versus untreated samples within each condition. For each trained condition (SYKi, MDP, BG), treatment-induced changes were compared with those observed in the untrained baseline by computing a delta-of-delta log2 fold change:

```text
Δlog2FC = log2FC(PMA-treated vs. untreated in trained condition) − log2FC(PMA-treated vs. untreated in untrained condition)
```

This approach identifies genes whose stimulation response is increased or decreased in trained conditions relative to the untrained baseline. Genes detected in both the trained and untrained comparisons were retained for ranking without additional significance filtering.

Ranked gene lists were analyzed using Gene Set Enrichment Analysis (GSEA) implemented in the clusterProfiler package against Blood Transcriptional Modules (BTM). P-values were adjusted using the Benjamini–Hochberg method.

Pathway enrichment results were generated for each cell type and training condition and summarized to compare condition-specific enrichment patterns across SYKi, MDP, BG, and the untrained baseline.

Filtered GSEA results were visualized as bubble plots using ggplot2. Pathways with adjusted p-values less than 0.05 were retained, and entries containing “TBA” in the pathway ID were excluded. Bubble plots were generated separately for MDP/BG, SYKi, and all conditions combined. In each plot, the x-axis represents training condition, the y-axis represents enriched pathway, bubble color represents normalized enrichment score (NES), and bubble size represents statistical significance as `-log10(p.adjust)`. Results were faceted by cell type using the `monaco_ann1` annotation, with cell types displayed in a manually defined order. Plots were exported as both PNG and SVG files for manuscript figure generation.

### ATAC-seq Associated Gene Pathway Analysis

Differential chromatin accessibility analysis revealed that statistically significant peaks passing multiple testing correction (FDR-adjusted threshold) were most robustly detected in β-glucan (BG)-trained CD14⁺ monocytes, whereas other training conditions and cell types showed fewer peaks meeting stringent significance thresholds. To enable comparative pathway-level interpretation across all training conditions, a peak-to-gene linkage–based enrichment approach was applied using correlated accessibility changes, allowing inclusion of biologically relevant regulatory signals that may not reach peak-level significance in all conditions. This approach is particularly suited to single-cell ATAC-seq data, where peak-level statistical power can be limited due to data sparsity.

Pathway enrichment analysis was performed on genes linked to differentially accessible ATAC-seq peaks using the clusterProfiler and ReactomePA packages. Peak-to-gene associations were derived from ArchR-generated correlation-based linkages between chromatin accessibility and gene activity scores. Peak–gene pairs were filtered using thresholds on differential accessibility (log2 fold change) and peak–gene correlation strength to retain high-confidence regulatory associations. Gene symbols were converted to ENTREZ identifiers using the org.Hs.eg.db annotation package.

Reactome and Gene Ontology (GO) enrichment analyses were performed on the resulting gene sets. Statistical significance was assessed using Benjamini–Hochberg correction, and pathways with adjusted p-values below the significance threshold were retained. Enrichment results were visualized using bar plots (Reactome, top pathways ranked by significance) and dot plots (GO, top pathways ranked by enrichment score). Analyses were conducted separately for each training condition (BG, MDP, SYKi) and cell type (CD14⁺ monocytes and CD4 TCM cells), focusing on positively correlated peak–gene associations indicative of activating regulatory relationships. Combined multi-panel visualizations were generated to compare shared and condition-specific pathway enrichments across training conditions and cell types.

### ATAC-seq Peak-to-Gene and Locus-Level Accessibility Analysis of ZNF Transcription Factors

To investigate regulatory changes associated with C2H2 zinc finger (ZNF) transcription factors, peak-to-gene linkage analysis was performed using ArchR in β-glucan (BG)-trained versus untrained conditions, where differential accessibility signals were most robust. Peak–gene links were generated using ArchR functions (`addPeak2GeneLinks`, `getPeak2GeneLinks`) and filtered to retain high-confidence associations based on correlation (>0.2), linkage false discovery rate (FDR < 0.05), peak-level differential accessibility (FDR < 0.05), and absolute log2 fold change (>0.5).

ZNF-associated peak–gene links were identified by selecting genes with names matching the “ZNF” prefix. To reduce redundancy arising from multiple gene assignments per peak, peak-to-gene links were ranked by linkage significance, and the most confident gene association was retained for each unique peak.

To distinguish direct locus accessibility from distal regulatory associations, genomic coordinates for ZNF genes were retrieved using biomaRt, extended by ±2 kb to capture promoter-proximal regions, and overlapped with ATAC-seq peaks using GenomicRanges. Peaks overlapping these regions were classified as locus-proximal, while non-overlapping peaks were considered distal regulatory elements.

Differential accessibility at ZNF loci was defined based on log2 fold change and FDR thresholds. This approach enables separation of direct chromatin accessibility changes at ZNF gene loci from distal regulatory interactions, providing insight into potential regulatory mechanisms associated with trained immunity–induced chromatin remodeling.

### Transcription Factor Motif Activity and Family Analysis

Transcription factor (TF) motif activity was quantified using chromVAR deviation scores implemented in ArchR. Motif deviation scores were computed from the ArchR MotifMatrix and represent GC-corrected, depth-normalized chromatin accessibility at peaks containing each motif. Differences in motif activity between conditions were assessed within each cell type using Wilcoxon rank-sum tests, and p-values were adjusted using the Benjamini–Hochberg procedure (FDR < 0.05).

Significant TF motifs identified from this analysis were mapped to transcription factor families using the CisBP database. For family-level analyses, TFs were aggregated irrespective of direction of change, such that both motifs with increased and decreased activity were included. This approach was used to assess overall TF family representation associated with differential chromatin accessibility, rather than direction-specific regulatory effects.

To evaluate overlap of TFs across conditions and cell types, intersection analyses were performed using precomputed TF sets and visualized with UpSet plots. TF family composition within each intersection was further summarized and visualized using bar plots and stacked bar charts to highlight dominant TF families across shared and condition-specific TF sets.

### Feature Visualization of Gene Expression and Module Scores

Gene expression and module score distributions were visualized across training conditions using Seurat-based FeaturePlot and violin plot representations. Cells were grouped by training condition and treatment status (untrained, SYKi, MDP, BG; ± PMA) based on sequencing lane metadata.

For module score visualization, precomputed module scores (e.g., `ModuleScore*`) stored in the Seurat metadata were plotted across conditions. For gene-level analysis, normalized expression values were extracted from the RNA assay and visualized using FeaturePlot projections on a shared low-dimensional embedding. If not already present, a two-dimensional embedding (UMAP or equivalent) was generated or reused for visualization.

For each feature, condition-specific FeaturePlots were generated using a consistent color scale, with expression values truncated at upper quantiles to improve visualization of dynamic range. Cell-type labels were overlaid using cluster centroids to aid interpretation.

Complementary violin plots were generated to summarize feature distributions across conditions. These plots include kernel density representations, overlaid boxplots, and median value annotations. For gene expression analyses, violin plots were optionally restricted to non-zero expression values to better highlight active expression patterns. Reference median values from the untrained condition were indicated using dashed lines to facilitate comparison across training conditions.

### Genome-wide Circos Visualization of ATAC Accessibility and HIV Integration

Genome-wide chromatin accessibility and HIV integration events were visualized using circos plots generated with the circlize package. Genomic bins containing accessibility counts for training conditions (BG, MDP, SYKi) were plotted as radial bar tracks across chromosomes. For selected analyses, HIV integration sites were overlaid as genomic points based on insertion coordinates, and differential ATAC-seq peaks were plotted as vertical line tracks filtered by log2 fold change (Log2FC > 1.5), correlation (> 0.3), and statistical significance (FDR ≤ 0.05).

Chromosome ordering was standardized to canonical human chromosomes (chr1–chr22, X, Y, M), and all tracks were aligned to a shared genomic coordinate system. Separate circos plots were generated for different cell types (e.g., CD4 TCM, CD14 monocytes), enabling comparison of genome-wide accessibility patterns across conditions. All plots were exported as PDF and SVG files for manuscript visualization.

### Genome Browser Visualization

HIV genome browser-style visualization was generated by mapping ATAC-seq peak accessibility signals and HIV gene annotations onto the HIV reference genome. Peak coordinates and scores were converted to genomic ranges, then stratified by training condition and PMA treatment status. For each condition, overlapping peak scores were summarized across genomic positions to generate ATAC-seq accessibility tracks. HIV gene annotations from the GFF3 file were plotted as a reference track, and condition-specific ATAC tracks were aligned to the same genomic coordinate range. Final multi-track genome browser plots were exported as PDF and SVG for manuscript visualization. The results are descriptive and used for visualization, not statistical testing.

### Chromatin Peak–HIV Integration Overlap Framework

Chromatin accessibility profiles were processed using ArchR, and peak–gene regulatory links were computed using `addPeak2GeneLinks` based on IterativeLSI embeddings. These links represent correlations between chromatin accessibility at ATAC peaks and gene activity, and only associations passing an FDR ≤ 0.05 threshold were retained.

To characterize accessibility differences between training conditions (BG, MDP, SYKi) and untrained cells, differential peak analysis was performed using `getMarkerFeatures` with a Wilcoxon test, controlling for TSSEnrichment and fragment count bias. To enable comprehensive assessment of potential HIV integration within accessible chromatin, peaks were not restricted to those passing differential accessibility FDR thresholds; instead, a minimal log2 fold change cutoff was applied to exclude low-signal regions while retaining a broad set of accessible loci.

For each training condition, peak coordinates, accessibility metrics, and peak–gene link information were combined to generate condition-specific peak-to-gene tables. HIV integration coordinates were intersected with these peak regions, with overlap defined as insertion sites falling within peak boundaries (peak_start ≤ insertion ≤ peak_end).

This framework enabled quantification of HIV integration events occurring within accessible chromatin regions linked to putative regulatory elements across conditions. Peak–gene associations were further stratified by correlation strength (positive: > 0.3; negative: < −0.3), and corresponding gene sets were analyzed using Reactome and Gene Ontology enrichment via clusterProfiler and ReactomePA.

Together, this approach integrates chromatin accessibility, regulatory linkage, and HIV integration site mapping to assess associations between accessible regulatory elements and patterns of viral integration across trained and untrained conditions.

### ZNF Peak-to-Gene and Local Locus Accessibility Analysis

Peak-to-gene (peak2gene) links derived from ATAC-seq data were filtered to identify peaks associated with zinc finger (ZNF) genes based on gene name matching. Peak-gene links were evaluated using correlation strength, link significance (FDR.y < 0.05), peak-level significance (FDR.x < 0.05 when available), and accessibility changes (|Log2FC| > 0.5).

To distinguish distal regulatory associations from physically proximal chromatin accessibility, genomic coordinates for ZNF genes were retrieved from Ensembl using biomaRt. Peaks were then intersected with ZNF gene bodies extended by ±2 kb using GenomicRanges, defining a set of locally overlapping peaks.

For each unique peak, the most significant peak-to-gene link was retained based on link FDR, correlation magnitude, and accessibility change. Peaks were categorized as either “local” (overlapping ZNF loci ±2 kb) or “non-local” (distal but linked). Summary statistics, gene-level aggregation, and visualization (bar plots, dot plots, and filter progression plots) were generated to compare peak2gene-linked versus physically local chromatin accessibility patterns.

```
```
