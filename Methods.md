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

Raw 10x Genomics scRNA-seq data were imported into Seurat objects per sequencing lane [(1,2)](#references), with sample identities assigned using Demuxalot [(6)](#references). Putative doublets were detected using DoubletFinder [(3)](#references) and removed prior to downstream processing. Lane-level objects were normalized, variable features identified, and merged into a single Seurat object for integration [(1,2)](#references). Batch correction and integration were performed with Seurat's CCA, RPCA [(1,2)](#references) and Harmony [(4)](#references) to support visualization and cell type annotation. Cell type annotation was performed using reference-based mapping with the Monaco immune reference dataset [(5)](#references) via Seurat's `FindTransferAnchors` and `TransferData` functions, followed by manual curation based on canonical marker gene expression to refine major immune populations and memory subsets. The final annotated object was used for all downstream analyses.

---

## Processing of scATAC-seq data

Single-cell ATAC-seq data were processed using the ArchR framework [(7)](#references). Raw fragment files were imported into ArchR projects, and quality control filtering was applied to remove low-quality nuclei based on transcription start site (TSS) enrichment and fragment counts. Putative doublets were identified and removed using ArchR's doublet detection algorithm [(7)](#references).

Dimensionality reduction was performed using iterative latent semantic indexing (LSI), followed by graph-based clustering and UMAP visualization [(7)](#references). Peaks were called using MACS2 within ArchR on pseudobulk profiles generated from clustered cells [(8)](#references).

To assign cell identities, gene activity scores were computed from chromatin accessibility profiles using ArchR [(7)](#references). Cell type labels were transferred from the matched scRNA-seq dataset to scATAC-seq cells using ArchR's integration framework [(7)](#references), enabling consistent annotation across modalities. Transferred labels were validated using accessibility at canonical marker genes.

Peak-to-gene linkages were computed using ArchR's `addPeak2GeneLinks` function based on IterativeLSI embeddings, retaining only associations with FDR ≤ 0.05 [(7)](#references). Motif enrichment analysis was performed using ArchR's motif annotation framework [(7)](#references), and transcription factor activity was inferred using chromVAR deviation scores [(9)](#references).

---

## Pathway Enrichment Analysis (scRNA-seq)

Single-cell RNA-seq data were aggregated into pseudobulk expression profiles using Seurat's `AggregateExpression()` function [(1,2)](#references), grouping by treatment status, cell type (`monaco_ann1` annotation), subject (`sCALL`), and sequencing lane to preserve biological and technical structure. Odd-numbered lanes (1, 3, 5, 7) were designated as untreated; even-numbered lanes (2, 4, 6, 8) as PMA-treated. Training conditions were defined by lane groups: Untrained (lanes 1, 2), SYKi (3, 4), MDP (5, 6), and BG (7, 8).

Differential expression analysis was performed on pseudobulk profiles using DESeq2 via Seurat's `FindMarkers` function (`test.use = "DESeq2"`), with a minimum detection threshold of `min.pct = 0.25`, comparing PMA-treated versus untreated samples within each condition. For each trained condition (SYKi, MDP, BG), treatment-induced changes were compared with those observed in the untrained baseline by computing a delta-of-delta log2 fold change:

$$\Delta\log_2\text{FC} = \log_2\text{FC}(\text{PMA-treated vs. untreated in trained}) - \log_2\text{FC}(\text{PMA-treated vs. untreated in untrained})$$

This approach identifies genes whose stimulation response is increased or decreased in trained conditions relative to the untrained baseline. Only genes detected in both the trained and untrained comparisons were retained for ranking; no additional significance filtering was applied prior to ranking to preserve the full gene list for GSEA.

Ranked gene lists were analyzed using Gene Set Enrichment Analysis (GSEA) implemented in the clusterProfiler package [(10)](#references) against Blood Transcriptional Modules (BTM), with gene set size limits of `minGSSize = 10` and `maxGSSize = 500`, an exponent of 1, and a permissive `pvalueCutoff = 2` to retain all results prior to filtering. P-values were adjusted using the Benjamini–Hochberg method.

Pathway enrichment results were generated for each cell type and training condition. Filtered GSEA results (adjusted p-value < 0.05; pathways with "TBA" in the pathway ID excluded) were visualized as bubble plots using ggplot2 [(15)](#references), in which the x-axis represents training condition, the y-axis represents enriched pathway, bubble color represents normalized enrichment score (NES) on a diverging scale (blue–white–red, midpoint = 0), and bubble size represents statistical significance as −log10(p.adjust). Separate plots were generated for MDP/BG, SYKi, and all conditions combined. Results were faceted by cell type using the `monaco_ann1` annotation. Plots were exported as PNG and SVG files.

---

## ATAC-seq Associated Gene Pathway Analysis

Differential chromatin accessibility analysis revealed that statistically significant peaks passing multiple testing correction (FDR-adjusted threshold) were most robustly detected in β-glucan (BG)-trained CD14⁺ monocytes, whereas other training conditions and cell types showed fewer peaks meeting stringent significance thresholds. To enable comparative pathway-level interpretation across all training conditions, a peak-to-gene linkage–based enrichment approach was applied using correlated accessibility changes, allowing inclusion of biologically relevant regulatory signals that may not reach peak-level significance in all conditions. This approach is particularly suited to single-cell ATAC-seq data, where peak-level statistical power can be limited due to data sparsity.

Pathway enrichment analysis was performed on genes linked to differentially accessible ATAC-seq peaks using the clusterProfiler [(10)](#references) and ReactomePA [(11)](#references) packages. Peak-to-gene associations were derived from ArchR-generated correlation-based linkages [(7)](#references) and filtered using the following thresholds: Log2FC > 1.5, peak–gene correlation > 0.3 (positive) or < −0.3 (negative), and link FDR ≤ 0.05. Peak-level FDR filtering (FDR.x) was not applied to maximize sensitivity across conditions. Gene symbols were converted to ENTREZ identifiers using the `org.Hs.eg.db` annotation package [(12)](#references).

Reactome pathway enrichment (`enrichPathway`, `pvalueCutoff = 0.05`) and Gene Ontology enrichment across all GO domains (`enrichGO`, `ont = "ALL"`, `pAdjustMethod = "BH"`, `pvalueCutoff = 0.05`, `qvalueCutoff = 0.2`) were performed on the resulting gene sets. Enrichment results were visualized using bar plots (Reactome, top 20 pathways) and dot plots (GO, top 20 pathways) using ggplot2 [(15)](#references). Analyses were conducted separately for each training condition (BG, MDP, SYKi) and cell type (CD14⁺ monocytes and CD4 TCM cells), focusing primarily on positively correlated peak–gene associations indicative of activating regulatory relationships. Combined multi-panel visualizations were exported as PDF and SVG.

---

## ATAC-seq Peak-to-Gene and Locus-Level Accessibility Analysis of ZNF Transcription Factors

To investigate regulatory changes associated with C2H2 zinc finger (ZNF) transcription factors, peak-to-gene linkage analysis was performed using ArchR [(7)](#references) in β-glucan (BG)-trained versus untrained CD14⁺ monocytes, where differential accessibility signals were most robust. Peak–gene links were generated using ArchR's `addPeak2GeneLinks` and `getPeak2GeneLinks` functions and filtered to retain high-confidence associations based on the following criteria: absolute correlation > 0.2, link FDR (FDR.y) < 0.05, peak-level FDR (FDR.x) < 0.05 (where available), and absolute log2 fold change > 0.5.

ZNF-associated peak–gene links were identified by selecting genes with names matching the `ZNF` prefix. To reduce redundancy arising from multiple gene assignments per peak, peak-to-gene links were ranked by linkage significance (FDR.y), then by absolute correlation and absolute log2 fold change, and the most confident gene association was retained for each unique peak using this priority ordering.

To distinguish direct locus accessibility from distal regulatory associations, genomic coordinates for ZNF genes were retrieved from Ensembl using biomaRt [(13)](#references) and extended by ±2 kb to capture promoter-proximal regions. These windows were overlapped with ATAC-seq peaks using GenomicRanges [(14)](#references). Peaks overlapping these regions were classified as locus-proximal; non-overlapping peaks were considered distal regulatory elements. Summary statistics, gene-level aggregation, and visualizations (bar plots, dot plots, and filter-progression plots) were generated to compare peak2gene-linked versus physically local chromatin accessibility patterns and exported as PDF and SVG.

---

## Transcription Factor Motif Activity and Family Analysis

Transcription factor (TF) motif activity was quantified using chromVAR deviation scores [(9)](#references) implemented in ArchR [(7)](#references). Motif deviation scores were computed from the ArchR `MotifMatrix` and represent GC-corrected, depth-normalized chromatin accessibility at peaks containing each motif. Differences in motif activity between conditions were assessed within each cell type using Wilcoxon rank-sum tests applied to pseudobulk-aggregated deviation scores, and p-values were adjusted using the Benjamini–Hochberg procedure (FDR < 0.05).

Significant TF motifs were mapped to transcription factor families using the CisBP database [(19)](#references). For family-level analyses, TFs were aggregated irrespective of direction of change, such that both motifs with increased and decreased activity were included, to assess overall TF family representation associated with differential chromatin accessibility.

To evaluate overlap of TFs across conditions and cell types, intersection analyses were performed using precomputed TF sets and visualized with UpSet plots using ComplexHeatmap [(17)](#references). TF family composition within each intersection was further summarized and visualized using bar plots and stacked bar charts generated with ggplot2 [(15)](#references). All outputs were exported as PDF and SVG.

---

## Feature Visualization of Gene Expression and Module Scores

Gene expression and module score distributions were visualized across training conditions using Seurat-based `FeaturePlot` and violin plot representations [(1,2)](#references). Cells were grouped by training condition and treatment status (untrained, SYKi, MDP, BG; ± PMA) based on sequencing lane metadata. Conditions were color-coded as: untrained (light blue), SYKi (blue), MDP (purple), BG (pink).

For module score visualization, precomputed module scores stored in the Seurat metadata were plotted across conditions. For gene-level analysis, normalized expression values were extracted from the RNA assay and projected onto a shared UMAP embedding, with expression values truncated at the 95th percentile to improve visualization of dynamic range. Cell-type labels were overlaid using cluster centroids.

Complementary violin plots were generated to summarize feature distributions across conditions, including kernel density representations, overlaid boxplots, and median value annotations. For gene expression analyses, violin plots were restricted to non-zero expression values to highlight active expression patterns. Reference median values from the untrained condition were indicated using dashed lines. All plots were exported as PDF and SVG files.

---

## Genome-wide Circos Visualization of ATAC Accessibility and HIV Integration

Genome-wide chromatin accessibility and HIV integration events were visualized using circos plots generated with the circlize package [(16)](#references). Genomic bins containing accessibility counts for training conditions (BG, MDP, SYKi) were plotted as radial bar tracks across chromosomes, color-coded as BG (red), SYKi (blue), and MDP (green). For analyses incorporating HIV integration sites, insertion coordinates were plotted as point tracks using midpoint positions within each bin, and differential ATAC-seq peaks were plotted as vertical line tracks filtered by Log2FC > 1.5, correlation > 0.3, and FDR ≤ 0.05.

Chromosome ordering was standardized to canonical human chromosomes (chr1–chr22, X, Y, M), and all tracks were aligned to a shared genomic coordinate system. Separate circos plots were generated for CD4 TCM and CD14⁺ monocytes. All plots were exported as PDF and SVG. These analyses are descriptive and were used for visualization only; no statistical testing was performed.

---

## Genome Browser Visualization

HIV genome browser-style visualization was generated using ggbio by mapping ATAC-seq peak accessibility signals and HIV gene annotations onto the HIV reference genome (HXB2; GFF3 annotation imported via rtracklayer). Peak coordinates and scores were converted to GRanges objects and stratified by training condition and PMA treatment status. For each condition, overlapping peak scores were summed across genomic positions to generate ATAC-seq accessibility tracks, which were aligned to a shared genomic coordinate range alongside the HIV gene annotation reference track. Final multi-track genome browser plots were exported as PDF and SVG. These results are descriptive and used for visualization; no statistical testing was performed.

---

## Chromatin Peak–HIV Integration Overlap Framework

Chromatin accessibility profiles were processed using ArchR [(7)](#references), and peak–gene regulatory links were computed using `addPeak2GeneLinks` based on IterativeLSI embeddings, retaining only associations with FDR ≤ 0.05.

To characterize accessibility differences between training conditions (BG, MDP, SYKi) and untrained cells, differential peak analysis was performed using `getMarkerFeatures` with a Wilcoxon test, controlling for `TSSEnrichment` and `log10(nFrags)` as bias variables. To enable comprehensive assessment of potential HIV integration within accessible chromatin, peaks were not restricted to those passing differential accessibility FDR thresholds; instead, a minimal Log2FC cutoff of > 0 was applied to exclude low-signal regions while retaining a broad set of accessible loci.

HIV integration site coordinates (derived from published integration site databases; see Data Availability) were intersected with peak regions, with overlap defined as insertion sites falling within peak boundaries (`peak_start ≤ insertion ≤ peak_end`). For each training condition, peak coordinates, accessibility metrics, and peak–gene link information were combined into condition-specific peak-to-gene tables.

Peak–gene associations were stratified by correlation strength (positive: > 0.3; negative: < −0.3), and corresponding gene sets were analyzed using Reactome and Gene Ontology enrichment via clusterProfiler [(10)](#references) and ReactomePA [(11)](#references). This framework integrates chromatin accessibility, regulatory linkage, and HIV integration site mapping to assess associations between accessible regulatory elements and patterns of viral integration across trained and untrained conditions.

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
19. Weirauch, M.T., et al., Determination and inference of eukaryotic transcription factor sequence specificity. *Cell*, 2014. 157(6): p. 1389–1401.
20. Wu, T., et al., clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*, 2021. 2(3): p. 100141.
