# Research Analysis Scripts

This repository contains a collection of R scripts for multi-omics data analysis and visualization. This serves as the starting point for generating a reproducible and trackable set of scripts and data files for experimental research.

## Scripts Overview

All analysis scripts are located in the `SCRIPTS/` directory, with outputs generated in the `MS_V1/` directory.

### Analysis Scripts

- **V1_pathways.R** - Generates pathway enrichment bubble plots for multiple cell types and conditions, saved as PNG and SVG formats in `MS_V1/pathway/`

- **V1_genome_browser.R** - Creates genome browser tracks combining ATAC-seq and scRNA-seq data across HIV genome positions, outputs PDF and SVG to `MS_V1/Genome_Browser/`

- **V1_TF_family.R** - Analyzes transcription factor families from intersection data, produces stacked bar charts, UpSet plots, and percentage-based visualizations in `MS_V1/TF/`

- **V1_Feature.R** - Generates Seurat-style feature plots and violin plots for gene expression and module scores across experimental conditions, saves to `MS_V1/feature_plot/`

- **V1_Circos.R** - Creates circular genome plots (circos plots) showing genomic bin data, HIV insertions, and ATAC peaks across chromosomes in `MS_V1/CIRCOS/`

- **V1_ATAC_Pathways.R** - Performs pathway enrichment analysis on ATAC-seq peak-associated genes, generating Reactome and GO enrichment plots in `MS_V1/ATAC_pathways/`

## Output Structure

```
MS_V1/
├── pathway/          # Pathway enrichment bubble plots
├── Genome_Browser/   # Genome browser visualization tracks  
├── TF/              # Transcription factor family analysis
├── feature_plot/    # Feature and violin plots
├── CIRCOS/          # Circular genome plots
└── ATAC_pathways/   # ATAC-seq pathway enrichment
```

## File Formats

All visualizations are generated in both PDF and SVG formats to ensure compatibility across different publication and presentation platforms.

File Formats
All visualizations are generated in both PDF and SVG formats to ensure compatibility across different publication and presentation platforms.
SVG Advantages for Publication
The SVG format provides significant advantages for figure preparation:

Easy text editing - Open SVG files directly in Adobe Illustrator to modify font sizes, labels, and annotations without quality loss
Scalable graphics - Vector format allows infinite scaling for different publication requirements (posters, papers, presentations)
Layer manipulation - Individual plot elements can be selected and modified independently (colors, line weights, transparency)
Publication-ready - Direct compatibility with journal submission systems and professional design software

## Reproducibility

This repository is designed to facilitate reproducible research by providing:
- Standardized output directory structure
- Consistent file naming conventions
- Both vector (SVG/PDF) and raster output formats
- Modular script organization for independent analysis components
