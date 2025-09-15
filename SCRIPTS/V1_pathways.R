# ---- Libraries ----
library(ggplot2)
library(dplyr)
library(stringr)
# library(svglite)  # <- only needed if svglite isn't already installed

# ---- Paths & Output Dir ----
combined_results_path <- file.path("./CURRENT_DATA/", "Filtered_GSEA_Results_v3.csv")
message("\nCombined GSEA results saved to ", combined_results_path, "\n")

out_dir <- "./MS_V1/pathway"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- Read ----
df <- read.csv(combined_results_path, stringsAsFactors = FALSE)

# ---- Helpers ----
custom_order <- c(
  "Monocytes", "T cells", "CD4+ T cells", "CD8+ T cells",
  "CD4 TCM", "CD8 TCM", "NK cells", "Dendritic cells", "B cells"
)

make_plot <- function(dat) {
  dat$CellType <- factor(dat$CellType, levels = custom_order)
  
  ggplot(dat, aes(x = Condition, y = Description, color = NES, size = -log10(p.adjust))) +
    geom_point(alpha = 0.8) +
    scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
    scale_size_continuous(range = c(4, 26)) +
    labs(
      x = "Condition", y = "Pathway",
      color = "NES", size = "-log10(p.adjust)",
      title = "Pathway Enrichment Analysis"
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16),
      axis.text.y = element_text(size = 18),
      strip.text  = element_text(size = 16, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 14),
      panel.spacing = unit(3, "lines")
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 3) +
    coord_cartesian(clip = "off")
}

# save_both <- function(plot, basename, width = 24, height = 34, dpi = 300) {
#   png_path <- file.path(out_dir, paste0(basename, ".png"))
#   svg_path <- file.path(out_dir, paste0(basename, ".svg"))
#   ggsave(png_path, plot, width = width, height = height, dpi = dpi)
#   ggsave(svg_path, plot, width = width, height = height, device = "svg")
#   message("Saved: ", basename, " (.png & .svg) -> ", out_dir)
# }


save_both <- function(plot, basename, width = 24, height = 34, dpi = 300,
                      units = "in", out_dir = getwd()) {
  # ensure output dir exists
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  svg_path <- file.path(out_dir, paste0(basename, ".svg"))
  
  # force ggplot2's ggsave and a single, valid units string
  ggplot2::ggsave(filename = png_path, plot = plot,
                  width = width, height = height, units = units, dpi = dpi)
  
  # For SVG, either rely on extension or set device = "svg" (svglite recommended)
  ggplot2::ggsave(filename = svg_path, plot = plot,
                  width = width, height = height, units = units, device = "svg")
  
  message("Saved: ", basename, " (.png & .svg) -> ", out_dir)
}


# =========================
# 1) MDP + BG
# =========================
sig_ids_mdp_bg <- df %>%
  filter(Condition %in% c("MDP", "BG")) %>%
  filter(p.adjust < 0.05) %>%
  filter(!str_detect(ID, "TBA")) %>%
  distinct(ID) %>%
  pull(ID)

df_mdp_bg <- df %>%
  filter(ID %in% sig_ids_mdp_bg, Condition %in% c("MDP", "BG"))

plot_mdp_bg <- make_plot(df_mdp_bg)
save_both(plot_mdp_bg, "v03_Pathway_Bubble_MDP_BG")

# =========================
# 2) SYKi only
# =========================
sig_ids_syki <- df %>%
  filter(Condition %in% c("SYKi")) %>%
  filter(p.adjust < 0.05) %>%
  filter(!str_detect(ID, "TBA")) %>%
  distinct(ID) %>%
  pull(ID)

df_syki <- df %>%
  filter(ID %in% sig_ids_syki, Condition %in% c("SYKi"))

plot_syki <- make_plot(df_syki)
save_both(plot_syki, "v03_Pathway_Bubble_SYKi")

# =========================
# 3) ALL conditions (p.adjust < 0.05, no TBA)
# =========================
sig_ids_all <- df %>%
  filter(p.adjust < 0.05) %>%
  filter(!str_detect(ID, "TBA")) %>%
  distinct(ID) %>%
  pull(ID)

df_all <- df %>%
  filter(ID %in% sig_ids_all)

plot_all <- make_plot(df_all)
save_both(plot_all, "v03_Pathway_Bubble_ALL")

message("All plots saved (PNG + SVG) to: ", out_dir, "\n")
