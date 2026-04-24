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
