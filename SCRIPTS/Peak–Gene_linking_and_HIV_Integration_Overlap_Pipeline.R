# Load required libraries
library(ArchR)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)

# Set working directory
setwd("/data/chi/TEMP/FOO/TEMP1/TEMP/L2_ALL/ALLv2")

# Load ArchR project
proj <- readRDS("archr_motif_deviations_clusters2/archr_motif_deviation_project.rds")
proj <- ArchR::addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

# Extract and update metadata with group information
metadata <- as.data.frame(proj@cellColData)
metadata <- metadata %>%
  mutate(Group = case_when(
    grepl("SYKi", rownames(metadata)) ~ "SYKi",
    grepl("MDP", rownames(metadata)) ~ "MDP",
    grepl("BG", rownames(metadata)) ~ "BG",
    TRUE ~ "Untrained" # Default to Untrained if no patterns match
  ))
proj@cellColData <- DataFrame(metadata)

# Define parameters
#cell_type <- "CD14 Mono"
cell_type <- "CD4 TCM"
training_sets <- c("SYKi", "MDP", "BG")
fdr_threshold <- 0.05
bg_threshold <- 2
#training_set<-"BG"

# Initialize final storage for results
all_combined_results <- list()

#table(proj@cellColData$Group)
#BG       MDP      SYKi Untrained 
#29430     30017     29726     30549 

# Loop over each training set
for (training_set in training_sets) {
  cat("Processing training set:", training_set, "\n")
  
  # Subset project for the specific cell type
  proj_celltype <- subsetCells(
    ArchRProj = proj,
    cellNames = proj$cellNames[proj$Clusters2 == cell_type]
  )
  
  # Perform differential peak analysis
  markers_training_vs_untrained <- getMarkerFeatures(
    ArchRProj = proj_celltype,
    useMatrix = "PeakMatrix",
    groupBy = "Group",
    useGroups = training_set,
    bgdGroups = "Untrained",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  # Extract Log2FC, FDR, and other data
  log2fc <- unlist(assay(markers_training_vs_untrained, "Log2FC"))
  fdr <- unlist(assay(markers_training_vs_untrained, "FDR"))
  row_data <- as.data.frame(rowData(markers_training_vs_untrained))
  assay_data <- as.data.frame(assay(markers_training_vs_untrained))
  
  # Combine extracted data
  combined_df <- cbind(row_data, assay_data)
  combined_df$Log2FC <- log2fc
  combined_df$FDR <- fdr
  
  # Filter significant peaks
  significant_df <- combined_df %>%
    # filter(FDR <= fdr_threshold & Log2FC >= bg_threshold) %>%
    # filter(Log2FC >= bg_threshold) %>%
    mutate(Significant_Peak_Name = paste0(seqnames, "_", start, "_", end))
  
  
  #  proj <- addPeak2GeneLinks(
  #    ArchRProj = proj
  #  )
  # Retrieve Peak2GeneLinks
  peak2gene_links <- getPeak2GeneLinks(proj)$Peak2GeneLinks
  # Convert Peak2GeneLinks to a data frame
  p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
  p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
  p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
  p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
  
  
  
  
  # Retrieve Peak2Gene links
  #peak2gene_links <- metadata(proj@peakSet)$Peak2GeneLinks
  #peak2gene_links$geneName <- mcols(metadata(peak2gene_links)$geneSet)$name[peak2gene_links$idxRNA]
  #peak2gene_links$peakName <- paste0(metadata(peak2gene_links), "_", start(peak2gene_links), "_", end(peak2gene_links))[peak2gene_links$idxATAC]
  # (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
  # Convert Peak2GeneLinks to data frame
  #p2geneDF <- as.data.frame(peak2gene_links)
  
  # Combine significant peaks with Peak2Gene links
  
  p2geneDF <- as.data.frame(p2geneDF)
  
  combined_results <- significant_df %>%
    dplyr::inner_join(p2geneDF, by = c("Significant_Peak_Name" = "peakName"))
  
  # Store results
  all_combined_results[[training_set]] <- combined_results
  
  # Save results to file
  write.csv(
    combined_results,
    file = paste0("v3_significant_peaks_", training_set, "_CD14_Mono.csv"),
    row.names = FALSE
  )
}

# Summary statistics across training sets
for (training_set in training_sets) {
  results <- all_combined_results[[training_set]]
  num_peaks <- results %>% distinct(idx) %>% nrow()
  num_genes <- results %>% distinct(geneName) %>% nrow()
  avg_genes_per_peak <- results %>%
    group_by(idxATAC) %>%
    summarise(num_genes = n_distinct(geneName)) %>%
    summarise(avg_genes = mean(num_genes)) %>%
    pull(avg_genes)
  
  cat("\nTraining Set:", training_set, "\n")
  cat("Number of unique peaks:", num_peaks, "\n")
  cat("Number of unique genes:", num_genes, "\n")
  cat("Average number of genes per peak:", avg_genes_per_peak, "\n")
}



##############


library(dplyr)
library(readr)
library(UpSetR)
library(ggplot2)

library(dplyr)
library(biomaRt)
library(readr)

# Step 1: Load human gene coordinates using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_coordinates <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "external_gene_name"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)

write_tsv(gene_coordinates, "gene_coordinates.tsv")

# Clean and rename columns for consistency
gene_coordinates <- gene_coordinates %>%
  dplyr::rename(
    Chromosome = chromosome_name,
    Gene_Start = start_position,
    Gene_End = end_position,
    Gene_Name = external_gene_name
  ) %>%
  filter(Chromosome %in% c(1:22, "X", "Y"))  # Filter only standard chromosomes

# Add "chr" prefix for consistency with your data
gene_coordinates$Chromosome <- paste0("chr", gene_coordinates$Chromosome)

# Step 2: Read mm2.csv (insertion hits)
mm2_hits <- read_csv("./mmc2.csv")  # Your insertion data

# Step 3: Perform overlap/join to find insertions in genes
hits_with_genes <- mm2_hits %>%
  inner_join(
    gene_coordinates,
    by = "Chromosome"  # Match chromosome
  ) %>%
  filter(start >= Gene_Start & start <= Gene_End)  # Check if insertion falls within the gene region

# Step 4: View results
hits_with_genes

gene_list<-hits_with_genes$Gene_Name %>% unique()



###############################
###############################
###############################
###############################
# Step 1: Read and process data across thresholds

#print(combined_results)
#Threshold CD4_TCM_BG_PercentHits CD4_TCM_SYK_PercentHits CD4_TCM_MDP_PercentHits CD14_Mono_BG_PercentHits CD14_Mono_SYK_PercentHits
#1       0.0               89.67033                90.03663                89.89011                 88.35165                  89.96337
#2       0.1               89.67033                90.03663                89.89011                 88.35165                  89.96337
#3       0.2               73.26007                76.11722                73.40659                 69.45055                  76.48352
#4       0.3               51.42857                52.08791                54.21245                 48.64469                  57.87546
#5       0.4               33.77289                34.06593                36.70330                 32.96703                  39.48718
#6       0.5               19.12088                20.36630                23.15018                 20.65934                  25.42125
#CD14_Mono_MDP_PercentHits
#1                  87.76557
#2                  87.76557
#3                  69.81685
#4                  46.44689
#5                  29.08425
#6                  16.11722


#BG=PEAKv3/v3_significant_peaks_BG_CD4_TCM.csv
#SYK=PEAKv3/v3_significant_peaks_SYKi_CD4_TCM.csv
#MDP=PEAKv3/v3_significant_peaks_MDP_CD4_TCM.csv

# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)

# File paths
file_paths <- list(
  BG = "PEAKv3/v3_significant_peaks_BG_CD4_TCM.csv",
  SYK = "PEAKv3/v3_significant_peaks_SYKi_CD4_TCM.csv",
  MDP = "PEAKv3/v3_significant_peaks_MDP_CD4_TCM.csv",
  BG2 = "PEAKv3/v3_significant_peaks_BG_CD14_Mono.csv",
  SYK2 = "PEAKv3/v3_significant_peaks_SYKi_CD14_Mono.csv",
  MDP2 = "PEAKv3/v3_significant_peaks_MDP_CD14_Mono.csv"
  
)

# Define thresholds
thresholds <- seq(0, 0.5, by = 0.1)

# Example gene list for comparison
# Initialize a list to store results for each file
results <- list()

# Process each file
for (name in names(file_paths)) {
  file_path <- file_paths[[name]]
  percent_hits <- data.frame()  # To store percent hits for this file
  
  # Loop through thresholds
  for (threshold in thresholds) {
    # Read and filter data
    filtered_data <- read_csv(file_path) %>%
      filter(!is.na(start)) %>%
      filter(Log2FC > 0.5) %>%
      filter(Correlation > threshold) %>%
      filter(FDR.y <= 0.05) %>%
      dplyr::select(chr = seqnames, start = start, end = end, value = Log2FC, FDR = FDR.y, geneName) %>%
      unique()
    
    # Calculate percent hits
    total_genes <- length(gene_list)
    intersection_hits <- intersect(gene_list, filtered_data$geneName)
    percent_hit <- (length(intersection_hits) / total_genes) * 100
    
    # Store results
    percent_hits <- rbind(percent_hits, data.frame(Threshold = threshold, PercentHits = percent_hit))
  }
  
  # Add to results list
  results[[name]] <- percent_hits
}

# Combine results into a single table
combined_results <- Reduce(function(x, y) merge(x, y, by = "Threshold", all = TRUE), results)
colnames(combined_results) <- c("Threshold", "CD4_TCM_BG_PercentHits", "CD4_TCM_SYK_PercentHits", "CD4_TCM_MDP_PercentHits", 
                                "CD14_Mono_BG_PercentHits", "CD14_Mono_SYK_PercentHits", "CD14_Mono_MDP_PercentHits")

# Display the table
print(combined_results)

# Optionally write to a CSV
write.csv(combined_results, "v3_percent_hits_comparison.csv", row.names = FALSE)






