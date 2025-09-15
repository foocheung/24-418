# Load necessary libraries
library(ggbio)
library(GenomicRanges)
library(rtracklayer)

# Load the HIV genome, peak data, and gene annotations
hiv_genome <- import("sequence.fasta", format = "fasta")
peak_data <- read.table("output.bed", header = TRUE)
#peak_data <- read.table("test.txt", header = TRUE)
gff_file <- import("./sequence.gff3", format = "gff")
gene_annotations <- subset(gff_file, type == "gene")

# Check and rename columns for clarity
colnames(peak_data) <- c("chromosome", "start", "end", "name", "score", "PMA_Activation", "Training_condition")
peak_data$start <- as.numeric(peak_data$start)
peak_data$end <- as.numeric(peak_data$end)
peak_data$score <- as.numeric(peak_data$score)

# Convert peak data to GRanges object
gr_peak_data <- GRanges(
  seqnames = Rle(peak_data$chromosome),
  ranges = IRanges(start = peak_data$start, end = peak_data$end),
  PMA_Activation = peak_data$PMA_Activation,
  Training_condition = peak_data$Training_condition,
  score = peak_data$score
)

# Function to create cumulative scores starting at 0 for each position
create_cumulative_scores <- function(subset_atac) {
  df <- as.data.frame(subset_atac)
  total_counts <- data.frame(start = integer(), total_count = numeric())
  
  # Loop through each position
  for (pos in seq(min(df$start), max(df$end))) {
    overlapping_rows <- df[df$start <= pos & df$end > pos, ]
    total_count <- sum(overlapping_rows$score)
    total_counts <- rbind(total_counts, data.frame(start = pos, total_count = total_count))
  }
  
  total_counts$cumulative_count <- cumsum(total_counts$total_count)
  
  return(total_counts)
}

# Subset the ATAC-seq data by both Training_condition and PMA_Activation (Treatment/Untreated)
atac_tracks <- list()
max_y_limit <- 0  # Variable to hold the maximum y-value across all tracks

for (group in unique(peak_data$Training_condition)) {
  for (treatment in unique(peak_data$PMA_Activation)) {
    subset_atac <- gr_peak_data[gr_peak_data$Training_condition == group & gr_peak_data$PMA_Activation == treatment]
    
    # Create cumulative scores
    cumulative_data <- create_cumulative_scores(subset_atac)
    
    # Update the max_y_limit based on total_count
    max_y_limit <- max(max_y_limit, max(cumulative_data$total_count))
    
    # Create a ggplot for each subset using geom_col for plotting
    atac_tracks[[paste(group, treatment)]] <- ggplot(data = cumulative_data) +
    #  geom_col(aes(x = start, y = total_count), fill = "steelblue", alpha = 0.7) +  # Using total_count for plotting
      geom_col(aes(x = start, y = total_count), fill = "black", alpha = 0.7) +  # Using total_count for plotting
      labs(title = paste("ATAC-seq:", group, "-", treatment),
           x = "Genomic Position",
           y = NULL) +  # Y-axis label removed
      theme_minimal() +
      xlim(0, max(cumulative_data$start) + 1000) +  # Adjust x limits based on data
      ylim(0, max_y_limit + 1) +  # Set the same y limits for all plots
      theme(axis.title.y = element_blank(),  # Remove y-axis title
            axis.text.y = element_text(size = 10),    # Y-axis text size
            axis.text.x = element_text(size = 10),    # X-axis text size
            panel.grid.major.y = element_line(color = "gray", size = 0.5),  # Show major gridlines for y-axis
            panel.grid.minor.y = element_blank(),      # Hide minor gridlines
            axis.ticks.y = element_line(size = 0.5))    # Show y-axis ticks
  }
}

# Prepare scRNA-seq data (from file "HIV_genes.info")
scrna_data <- read.table("./HIV_genes.info", header = TRUE, sep = "\t")
gene_positions <- data.frame(
  gene = c("HIV1gp1", "HIV1gp3", "HIV1gp5", "HIV1gp8"),
  start = c(336, 4587, 5377, 5771),
  end = c(4642, 5165, 7970, 8341)
)
merged_data <- merge(scrna_data, gene_positions, by.x = "Gene", by.y = "gene")

# Create GRanges object for scRNA-seq data
scRNAseq_gr <- GRanges(
  seqnames = Rle("NC_001802.1"),
  ranges = IRanges(start = merged_data$start, end = merged_data$end),
  gene = merged_data$Gene,
  type = "gene",
  expression = merged_data$Count,  # Expression values
  group = merged_data$Group,
  treatment = merged_data$Treatment,
  genome = merged_data$Genome
)


# Subset scRNA-seq data by Group and create tracks



# Load necessary libraries
library(ggbio)
library(ggplot2)

# Function to safely add tracks if they exist or return an empty track if missing
add_track_if_exists <- function(track_list, track_name, xlim_range = NULL) {
  if (track_name %in% names(track_list)) {
    return(track_list[[track_name]] + theme(legend.position = "none") + 
             (if (!is.null(xlim_range)) xlim(xlim_range[1], xlim_range[2])))
  } else {
    # Return an empty placeholder track with the same x-axis limits
    empty_track <- ggplot() +
      geom_blank() +  # An empty plot
      theme_void() +  # No axes, labels, etc.
      labs(title = paste(track_name))  # Title for the track
    if (!is.null(xlim_range)) {
      empty_track <- empty_track + xlim(xlim_range[1], xlim_range[2])
    }
    return(empty_track)
  }
}


scRNAseq_tracks <- list()

# Find the maximum number of stacked boxes in any track
max_stack <- max(table(scRNAseq_gr$group, scRNAseq_gr$treatment))

for (group in unique(merged_data$Group)) {
  for (treatment in unique(merged_data$Treatment)) {
    subset_scrna <- scRNAseq_gr[scRNAseq_gr$group == group & scRNAseq_gr$treatment == treatment]
    
    # Skip if the subset is empty (no ranges)
    if (length(subset_scrna) == 0) next
    
    df <- as.data.frame(subset_scrna)
    
    # Sort by start position
    df <- df[order(df$start), ]
    
    # Create a new column for y-axis stacking
    df$stack_y <- ave(df$start, df$seqnames, FUN = function(x) seq_along(x))
    
    # Create the plot with stacked boxes
    scRNAseq_tracks[[paste(group, treatment)]] <- ggplot(df) +
    #  geom_tile(aes(x = (start + end) / 2, y = stack_y, width = end - start, height = 0.8, fill = genome), 
    #            color = "black") +  # Add a border around the boxes
       geom_tile(aes(x = (start + end) / 2, y = stack_y, width = end - start, height = 0.8), 
                 color = "black") +  # Add a border around the boxes
      
         #   scale_fill_manual(values = c("HIV-1" = "black", "HIV-M" = "blue")) +  # Example fill colors
      labs(title = paste("scRNA-seq:", group, "-", treatment), y = NULL) +
      scale_y_continuous(labels = NULL) +
      theme_minimal() +
      ylim(0, max_stack + 1) +  # Set consistent y-axis limits based on the maximum number of stacks
      theme(
        axis.text.y = element_blank(),     # Hide y-axis text
        axis.ticks.y = element_blank(),    # Hide y-axis ticks
        axis.title.y = element_blank(),    # Hide y-axis title
        panel.grid.major.y = element_blank(),  # Remove major y-axis gridlines
        panel.grid.minor.y = element_blank(),  # Remove minor y-axis gridlines
        axis.line.y = element_blank()      # Remove y-axis line itself
      )
  }
}



# Define xlim_range if not already defined
xlim_range <- c(min(start(gene_annotations)), max(end(gene_annotations)))


# Determine x-axis limits based on the gene annotation data
# xlim_range <- range(c(start(gene_annotations), end(gene_annotations)))
# 
# # Add gene annotations for HIV genome with fixed x-axis limits
# gene_annotations$gene_id <- gene_annotations$Name  
# gene_annotation_track <- autoplot(gene_annotations) + 
#   geom_text(aes(x = (start(gene_annotations) + end(gene_annotations)) / 2, 
#                 y = 0, 
#                 label = gene_annotations$gene_id),  # Label genes
#             size = 3, fontface = "bold", vjust = 0) +  # Text style
#   labs(title = "HIV Genes with IDs") + 
#   theme_minimal() +
#   scale_x_continuous(limits = c(xlim_range[1], xlim_range[2]))  # Fix x-axis to the genomic range
library(ggplot2)
library(ggrepel)

# Add gene annotations for HIV genome with fixed x-axis limits
gene_annotations$gene_id <- gene_annotations$Name

# Plot using autoplot and ggrepel without lines
gene_annotation_track<-autoplot(gene_annotations) + 
  geom_text_repel(aes(
    x = (start(gene_annotations) + end(gene_annotations)) / 2,
    y = -0.5,
    label = gene_annotations$gene_id
  ), 
  size = 5, fontface = "bold",# vjust = 0.005, 
  nudge_x = 0.005,       # Nudge text vertically to avoid overlap
  direction = "x",     # Only repel vertically
  segment.color = NA,  # Disable line segments
  max.overlaps = 5    # Max overlaps to consider
  ) + 
  labs(title = "HIV Genes with IDs") + 
  theme_minimal() +
  scale_x_continuous(limits = c(xlim_range[1], xlim_range[2]))


# Create final multi-track plot with ATAC-seq and scRNA-seq data split by treatment
p <- tracks(
  Genes = gene_annotation_track + theme(legend.position = "none"),
#  "ATAC\nUntrained\nUntreated" = add_track_if_exists(atac_tracks, "untrained Untreated", xlim_range),
  "ATAC\nUntrained\nTreated" = add_track_if_exists(atac_tracks, "untrained Treated", xlim_range),
#  "ATAC\nSYKi\nUntreated" = add_track_if_exists(atac_tracks, "SYKi Untreated", xlim_range),
  "ATAC\nSYKi\nTreated" = add_track_if_exists(atac_tracks, "SYKi Treated", xlim_range),
#  "ATAC\nMDP\nUntreated" = add_track_if_exists(atac_tracks, "MDP Untreated", xlim_range),
  "ATAC\nMDP\nTreated" = add_track_if_exists(atac_tracks, "MDP Treated", xlim_range),
#  "ATAC\nBG\nUntreated" = add_track_if_exists(atac_tracks, "BG Untreated", xlim_range),
  "ATAC\nBG\nTreated" = add_track_if_exists(atac_tracks, "BG Treated", xlim_range),
#  "scRNA\nUntrained\nUntreated" = add_track_if_exists(scRNAseq_tracks, "untrained Untreated", xlim_range),
#  "scRNA\nUntrained\nTreated" = add_track_if_exists(scRNAseq_tracks, "untrained Treated", xlim_range),
#  "scRNA\nSYKi\nUntreated" = add_track_if_exists(scRNAseq_tracks, "SYKi Untreated", xlim_range),
##  "scRNA\nSYKi\nTreated" = add_track_if_exists(scRNAseq_tracks, "SYKi Treated", xlim_range),
##  "scRNA\nMDP\nUntreated" = add_track_if_exists(scRNAseq_tracks, "MDP Untreated", xlim_range),
#  "scRNA\nMDP\nTreated" = add_track_if_exists(scRNAseq_tracks, "MDP Treated", xlim_range),
#  "scRNA\nBG\nUntreated" = add_track_if_exists(scRNAseq_tracks, "BG Untreated", xlim_range),
#  "scRNA\nBG\nTreated" = add_track_if_exists(scRNAseq_tracks, "BG Treated", xlim_range),
  title = "",
  label.width = unit(4.5, "lines")  # Adjust the label width
)

out_dir <- file.path("MS_V1", "Genome_Browser")

# Display the final plot
#print(p)

# Save as PDF
#pdf("genome_browser_6.pdf", width = 12, height = 5)
p + theme_tracks_sunset() +  #theme(
  #axis.text.y = element_blank(),     # Hide y-axis text
  #axis.ticks.y = element_blank(),    # Hide y-axis ticks
  #axis.title.y = element_blank(),    # Hide y-axis title
  #panel.grid.major.y = element_blank(),  # Remove major y-axis gridlines
  #panel.grid.minor.y = element_blank(),  # Remove minor y-axis gridlines
  #axis.line.y = element_blank(),
  #legend.position = "none" # Remove y-axis line itself
  #  )  
  guides(fill = "none") 
 # +  # Remove legend for fill aesthetic
  # theme_minimal() 

#dev.off()

grDevices::pdf(file.path(out_dir, "genome_browser_6.pdf"), width = 12, height = 5)
print(p)
dev.off()

# SVG
svglite::svglite(file.path(out_dir, "genome_browser_6.svg"), width = 12, height = 5)
print(p)
dev.off()