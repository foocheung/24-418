# Methods

## Table of Contents
- [Processing scRNA-seq](#processing-scrna-seq)
- [Processing of scATAC-seq data](#processing-of-scatac-seq-data)
- [Pathway Enrichment Analysis (scRNA-seq)](#pathway-enrichment-analysis-scrna-seq)
- [ATAC-seq Associated Gene Pathway Analysis](#atac-seq-associated-gene-pathway-analysis)
- [ATAC-seq Peak-to-Gene and Locus-Level Accessibility Analysis of ZNF Transcription Factors](#atac-seq-peak-to-gene-and-locus-level-accessibility-analysis-of-znf-transcription-factors)
- [Transcription Factor Motif Activity and Family Analysis](#transcription-factor-motif-activity-and-family-analysis)
- [Feature Visualization of Gene Expression and Module Scores](#feature-visualization-of-gene-expression-and-module-scores)
- [Genome-wide Circos Visualization of ATAC Accessibility and HIV Integration](#genome-wide-circos-visualization-of-atac-accessibility-and-hiv-integration)
- [Genome Browser Visualization](#genome-browser-visualization)
- [Chromatin Peak–HIV Integration Overlap Framework](#chromatin-peakhiv-integration-overlap-framework)
- [References](#references)

---

## Processing scRNA-seq

Raw 10x Genomics scRNA-seq data were imported into Seurat objects per sequencing lane (1,2), with sample identities assigned using Demuxalot (6). Putative doublets were detected using DoubletFinder (3) and removed prior to downstream processing. Lane-level objects were normalized, variable features were identified, and data were merged into a single Seurat object (1,2).

Batch correction and integration were performed using Seurat’s canonical correlation analysis (CCA) framework (1,2), and the resulting integrated embeddings (integrated.cca) were used for UMAP visualization and downstream cell type annotation. Alternative integration methods, including Harmony (4), were evaluated and compared. CCA was selected for downstream analyses as it provided improved preservation of known cell type structure and more consistent alignment across treated and untreated conditions.

Cell type annotation was performed using reference-based mapping with the Monaco immune reference dataset (5) via Seurat's FindTransferAnchors and TransferData functions, followed by manual curation based on canonical marker gene expression to refine major immune populations and memory subsets.

---

## Processing of scATAC-seq data

Single-cell ATAC-seq data were processed using the ArchR framework (7). Raw fragment files were imported into ArchR projects, and quality control filtering was applied to remove low-quality nuclei based on transcription start site (TSS) enrichment and fragment counts. Putative doublets were identified and removed using ArchR's doublet detection algorithm (7).

Dimensionality reduction was performed using iterative latent semantic indexing (LSI), followed by graph-based clustering and UMAP visualization (7). Peaks were called using MACS2 on pseudobulk accessibility profiles generated from clustered cells within ArchR (8).

To assign cell identities, gene activity scores were computed from chromatin accessibility profiles using ArchR (7). Cell type labels were transferred from the matched scRNA-seq dataset generated from the same samples using ArchR's integration framework (7), enabling consistent annotation across modalities. Transferred labels were validated using accessibility at canonical marker genes.

Peak-to-gene linkages were computed using ArchR's addPeak2GeneLinks function based on IterativeLSI embeddings, retaining associations with linkage FDR ≤ 0.05 for downstream analyses (7).

Motif enrichment analysis was performed using ArchR’s motif annotation framework, and transcription factor activity was quantified using chromVAR deviation scores, which account for GC bias and sequencing depth (9).

---

## Pathway Enrichment Analysis (scRNA-seq)

Single-cell RNA-seq data were aggregated into pseudobulk expression profiles using Seurat's AggregateExpression() function (1,2), grouping by treatment status, cell type (monaco_ann1 annotation), subject (sCALL), and sequencing lane to preserve biological and technical structure.

Treatment status (PMA-treated vs untreated) and training conditions (Untrained, SYKi, MDP, BG) were defined based on experimental design and encoded in sequencing metadata.

Differential expression analysis was performed using a DESeq2-based approach implemented in Seurat (FindMarkers, test.use = "DESeq2") (18). This analysis was performed specifically to generate log2 fold change values for each gene, which were used to rank genes for downstream pathway enrichment analysis, rather than for formal differential expression inference.

For each trained condition (SYKi, MDP, BG), treatment-induced changes were compared with those observed in the untrained baseline by computing a delta-of-delta log2 fold change:

Δlog2FC = log2FC (PMA-treated vs untreated in trained) − log2FC (PMA-treated vs untreated in untrained)

This approach isolates training-specific transcriptional responses by controlling for baseline stimulation effects observed in untrained cells. Only genes detected in both the trained and untrained comparisons were retained for ranking; no additional significance filtering was applied prior to ranking to preserve the full gene list for downstream enrichment analysis.

Ranked gene lists were analyzed using Gene Set Enrichment Analysis (GSEA) implemented in the clusterProfiler package (10,20) against Blood Transcriptional Modules (BTM), with gene set size limits of minGSSize = 10 and maxGSSize = 500, an exponent of 1, and a permissive pvalueCutoff = 2 to retain all pathways prior to downstream filtering. P-values were adjusted using the Benjamini–Hochberg method.

Pathway enrichment results were generated for each cell type and training condition. Filtered GSEA results (adjusted p-value < 0.05; pathways with "TBA" in the pathway ID excluded) were visualized as bubble plots using ggplot2 (15), where the x-axis represents training condition, the y-axis represents enriched pathway, bubble color represents normalized enrichment score (NES), and bubble size represents statistical significance as −log10(p.adjust). Results were faceted by cell type using the monaco_ann1 annotation and exported for downstream visualization.

---

## ATAC-seq Associated Gene Pathway Analysis

Thresholds for log2 fold change, peak–gene correlation, and statistical significance were adapted across analyses based on the specific biological question and the need to balance stringency and signal retention in sparse single-cell ATAC-seq data.

Differential chromatin accessibility analysis revealed that statistically significant peaks passing multiple testing correction were most robustly detected in β-glucan (BG)-trained CD14⁺ monocytes, whereas other training conditions and cell types showed fewer peaks meeting stringent significance thresholds. To enable pathway-level comparisons across all training conditions, a peak-to-gene linkage–based enrichment approach was applied, leveraging correlated accessibility changes to capture biologically relevant regulatory signals that may not reach peak-level statistical significance in all conditions.

Pathway enrichment analysis was performed on genes linked to accessible chromatin regions using the clusterProfiler (10,20) and ReactomePA (11) packages. Peak-to-gene associations were derived from ArchR-generated correlation-based linkages (7) and filtered using the following criteria: absolute log2 fold change > 1.5, peak–gene correlation > 0.3 (positive) or < −0.3 (negative), and linkage FDR ≤ 0.05. Peak-level FDR filtering was not applied in order to retain sufficient signal across conditions with lower statistical power.

Gene symbols were converted to ENTREZ identifiers using the org.Hs.eg.db annotation package (12). Reactome pathway enrichment (enrichPathway, pvalueCutoff = 0.05) and Gene Ontology enrichment across all domains (enrichGO, ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2) were performed on the resulting gene sets.

Analyses were conducted separately for each training condition (BG, MDP, SYKi) and cell type (CD14⁺ monocytes and CD4 TCM cells), focusing primarily on positively correlated peak–gene associations indicative of activating regulatory relationships. Enrichment results were visualized using bar plots (Reactome, top pathways) and dot plots (GO) generated with ggplot2 (15), and exported for downstream interpretation.

---

## ATAC-seq Peak-to-Gene and Locus-Level Accessibility Analysis of ZNF Transcription Factors

To investigate regulatory changes associated with C2H2 zinc finger (ZNF) transcription factors, peak-to-gene linkage analysis was performed using ArchR (7) in β-glucan (BG)-trained versus untrained CD14⁺ monocytes, where differential accessibility signals were most robust.

Peak–gene links were generated using ArchR’s addPeak2GeneLinks and getPeak2GeneLinks functions and filtered to retain high-confidence associations based on the following criteria: absolute correlation > 0.2, peak-to-gene linkage FDR < 0.05, and absolute log2 fold change > 0.5, with peak-level FDR < 0.05 applied where available. These thresholds were selected to capture a broader set of regulatory associations relative to pathway-level analyses while maintaining biological relevance.

ZNF-associated peak–gene links were identified by selecting genes with names matching the “ZNF” prefix. Because multiple peak–gene associations can exist for a given genomic region, peak–gene rows were ranked by linkage significance (FDR), followed by absolute correlation and absolute log2 fold change, and deduplicated by genomic coordinates to retain a single highest-confidence peak–gene association per peak.

To distinguish local chromatin accessibility from distal regulatory interactions, genomic coordinates for ZNF genes were retrieved from Ensembl using biomaRt (13) and extended by ±2 kb to define promoter-proximal regions. These regions were overlapped with ATAC-seq peaks using GenomicRanges (14). Peaks overlapping these regions were classified as locus-proximal, while non-overlapping peaks were classified as distal regulatory elements.

Summary statistics and gene-level aggregations were generated to compare peak-to-gene linked accessibility with locus-proximal accessibility. Results were visualized using bar plots and dot plots and exported for downstream interpretation.

---

## Transcription Factor Motif Activity and Family Analysis

Transcription factor (TF) motif activity was quantified using chromVAR deviation scores (9) implemented in ArchR (7). Motif deviation scores were computed from the ArchR MotifMatrix and represent GC-corrected, depth-normalized chromatin accessibility at peaks containing each motif.

Differences in motif activity between conditions were assessed within each cell type using Wilcoxon rank-sum tests applied to single-cell deviation scores. P-values were adjusted using the Benjamini–Hochberg procedure, and statistical significance was defined as FDR < 0.05.

Significant TF motifs were mapped to transcription factor families using the CisBP database (19). For family-level analyses, TFs were aggregated irrespective of direction of change, such that both motifs with increased and decreased activity were included to assess overall TF family representation associated with differential chromatin accessibility.

To evaluate overlap of TFs across conditions and cell types, intersection analyses were performed using precomputed TF sets and visualized with UpSet plots using ComplexHeatmap (17). TF family composition within each intersection was further summarized and visualized using bar plots and stacked bar charts generated with ggplot2 (15). All outputs were exported for downstream interpretation.

---

## Feature Visualization of Gene Expression and Module Scores

Gene expression and module scores were visualized using Seurat-based FeaturePlot and violin plot functions (1,2). Cells were grouped by training condition and treatment status (Untrained, SYKi, MDP, BG; ± PMA) based on sequencing metadata.

For gene expression visualization, normalized expression values from the RNA assay were projected onto a shared UMAP embedding generated from the integrated dataset. Expression values were capped at the 95th percentile to improve visualization of dynamic range.

Module scores precomputed in the Seurat metadata were visualized across conditions using both UMAP-based FeaturePlots and violin plots. Violin plots included kernel density distributions, overlaid boxplots, and median value annotations to summarize expression patterns across conditions.

For gene-level analyses, violin plots were restricted to non-zero expression values to highlight active expression patterns. Reference median values from the untrained condition were indicated using dashed lines where applicable.

All plots were generated using Seurat and ggplot2 (15) and exported for downstream visualization.

---

## Genome-wide Circos Visualization of ATAC Accessibility and HIV Integration

Genome-wide chromatin accessibility and HIV integration events were visualized using circos plots generated with the circlize package (16). To summarize genome-wide patterns, gene-level genomic distributions were binned across the human genome using fixed-width intervals of 1 Mb generated with the tileGenome function based on hg38 genome coordinates.

Gene sets associated with differential chromatin accessibility were first derived from peak-to-gene linkage analyses, and genomic coordinates for these genes were obtained using Ensembl annotations via biomaRt. These gene coordinates were converted into genomic ranges, and overlaps between genes of interest and genomic bins were quantified using the countOverlaps function, producing bin-level counts for each training condition (BG, MDP, SYKi).

Binned gene counts were visualized as radial bar tracks in circos plots to represent genome-wide accessibility-associated gene distributions. HIV integration sites were plotted as point tracks based on genomic coordinates, with insertion positions mapped to bin midpoints for visualization.

Differential ATAC-seq peaks were incorporated as additional genomic tracks and visualized as vertical line features scaled by log2 fold change, filtered using thresholds of log2 fold change > 1.5, peak–gene correlation > 0.3, and linkage FDR ≤ 0.05.

Although statistical enrichment of genomic bins was assessed using hypergeometric testing with multiple testing correction, circos visualizations display raw bin-level counts for descriptive purposes only. All genomic features were aligned to a common coordinate system across canonical human chromosomes (chr1–chr22, X, Y). Separate circos plots were generated for CD4 TCM cells and CD14⁺ monocytes.

---

## Genome Browser Visualization

Genome browser visualizations were generated to examine chromatin accessibility and peak-level features at specific genomic loci. Coverage tracks were produced from scATAC-seq fragment data using ArchR (7), with accessibility signal normalized for sequencing depth and visualized as continuous tracks across genomic regions of interest.

Differentially accessible peaks were overlaid as annotation tracks based on genomic coordinates and filtered using thresholds of log2 fold change > 1.5, peak–gene correlation > 0.3, and linkage FDR ≤ 0.05. Peak-to-gene linkages were incorporated to contextualize regulatory interactions at selected loci.

Genome browser tracks were visualized using standard genomic visualization tools, enabling comparison of accessibility patterns across training conditions (BG, MDP, SYKi) and treatment states. Representative loci were selected for visualization based on biological relevance and differential accessibility patterns.

---

## Chromatin Peak–HIV Integration Overlap Framework

To assess the relationship between chromatin accessibility and HIV integration, a gene-level overlap framework was applied integrating differential accessibility, peak-to-gene regulatory linkage, and genomic mapping of viral insertion sites.

Chromatin accessibility profiles were processed using ArchR (7), and peak–gene regulatory links were computed using addPeak2GeneLinks based on IterativeLSI embeddings. Differential peak analysis was performed using getMarkerFeatures with a Wilcoxon test, controlling for TSSEnrichment and log10(nFrags) as bias variables. Peak-level statistics, including log2 fold change and FDR, were extracted and combined with peak-to-gene linkage data by matching peak coordinates.

Peaks were filtered to retain accessible regulatory regions using thresholds of log2 fold change > 0.5 and peak-to-gene linkage FDR ≤ 0.05, with correlation thresholds (> 0 to > 0.5) varied to assess robustness of regulatory associations. These thresholds were selected to balance sensitivity and specificity and to evaluate the stability of gene-level overlap across a range of regulatory stringencies. No strict peak-level differential accessibility FDR cutoff was applied, allowing inclusion of a broader set of potentially relevant accessible regions.

HIV integration site coordinates (derived from published data) were mapped to gene bodies using genomic overlap with Ensembl gene annotations obtained via biomaRt (13), generating a reference gene set associated with viral insertion events.

For each training condition (BG, MDP, SYKi) and cell type, peak-associated gene sets derived from filtered peak-to-gene linkages were compared with HIV integration-associated gene sets. Overlap was quantified as the proportion of HIV-associated genes detected within peak-linked gene sets across varying correlation thresholds, enabling evaluation of the relationship between chromatin accessibility, regulatory linkage, and viral integration patterns.

---

## References

1. Stuart, T., et al., Comprehensive integration of single-cell data. *Cell*, 2019. 177(7): p. 1888–1902.
2. Hao, Y., et al., Integrated analysis of multimodal single-cell data. *Cell*, 2021. 184(13): p. 3573–3587.
3. McGinnis, C.S., Murrow, L.M., and Gartner, Z.J., DoubletFinder: Doublet detection in single-cell RNA sequencing data using artificial nearest neighbors. *Cell Systems*, 2019. 8(4): p. 329–337.
4. Korsunsky, I., et al., Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 2019. 16(12): p. 1289–1296.
5. Monaco, G., et al., RNA-Seq signatures normalized by mRNA abundance allow absolute deconvolution of human immune cell types. *Cell Reports*, 2019. 26(6): p. 1627–1640.
6. Rogozhnikov, A., et al., Demuxalot: scaled up genetic demultiplexing for single-cell sequencing. bioRxiv, 2021: p. 2021.05.
7. Granja, J.M., et al., ArchR: An integrative and scalable software package for single-cell chromatin accessibility analysis. *Nature Genetics*, 2021. 53(3): p. 403–411.
8. Zhang, Y., et al., Model-based analysis of ChIP-Seq (MACS). *Genome Biology*, 2008. 9(9): p. R137.
9. Schep, A.N., et al., chromVAR: Inferring transcription-factor-associated accessibility from single-cell epigenomic data. *Nature Methods*, 2017. 14(10): p. 975–978.
10. Yu, G., et al., clusterProfiler: An R package for comparing biological themes among gene clusters. *OMICS*, 2012. 16(5): p. 284–287.
11. Yu, G. and He, Q.Y., ReactomePA: An R/Bioconductor package for reactome pathway analysis and visualization. *Molecular BioSystems*, 2016. 12(2): p. 477–479.
12. Carlson, M., org.Hs.eg.db: Genome-wide annotation for Human. R package version 3.22.0, 2019.
13. Durinck, S., et al., BioMart and Bioconductor: A powerful link between biological databases and microarray data analysis. *Bioinformatics*, 2005. 21(16): p. 3439–3440.
14. Lawrence, M., et al., Software for computing and annotating genomic ranges. *PLoS Computational Biology*, 2013. 9(8): p. e1003118.
15. Wickham, H., ggplot2: Elegant graphics for data analysis. Springer-Verlag New York, 2016.
16. Gu, Z., et al., circlize implements and enhances circular visualization in R. *Bioinformatics*, 2014. 30(19): p. 2811–2812.
17. Lex, A., et al., UpSet: Visualization of intersecting sets. *IEEE Transactions on Visualization and Computer Graphics*, 2014. 20(12): p. 1983–1992.
18. Love, M.I., Huber, W., and Anders, S., Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 2014. 15(12): p. 550.
