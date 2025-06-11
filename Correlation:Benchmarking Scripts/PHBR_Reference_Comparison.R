# PHBR Comparison with Reference & Preprocessing Script
# --------------------------------------------
# Load necessary libraries
# --------------------------------------------
library(glue)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# --------------------------------------------
# Load core data
# --------------------------------------------
mut_matrix <- readRDS("../data/mut_matrix.rds")
data_per_cancer_type <- readRDS("../data/data_per_cancer_type.rds")
patient_info <- read.csv("../data/patient_info.csv", stringsAsFactors = FALSE)

# --------------------------------------------
# Set up PHBR file paths
# --------------------------------------------
phbr_folder <- ".../phbr_nextflow/results/phbr"
phbr_files <- list.files(phbr_folder, pattern = "_mhc1.csv$", full.names = TRUE)

# Log number of PHBR files
message(glue("Found {length(phbr_files)} PHBR files."))

# Extract patient IDs from filenames
patient_ids <- sub(".*filtered_mhc1_([^_]+)*", "\\1", basename(phbr_files))
num_patients <- length(unique(patient_ids))
cat("Number of unique patients:", num_patients, "\n")

# Example access to mutation data
m <- mut_matrix["TCGA-FI-A2CX-01", ]
m[m == 1]

# Example access to affinity data
data_per_cancer_type$BLCA$affinities["TCGA-XF-A9SZ",]["TP53_R156P"]

# --------------------------------------------
# Build lookup for cancer types
# --------------------------------------------
cancer_type_map <- setNames(patient_info$cancer_type, patient_info$sample_id)
message(glue("Loaded cancer types for {length(cancer_type_map)} samples."))

# --------------------------------------------
# Main PHBR processing loop
# --------------------------------------------

all_results <- list()

for (file in phbr_files) {
  
  # Extract sample ID from filename
  sample_id <- sub("phbr_(TCGA-[A-Z0-9-]+)_mhc1.csv", "\\1", basename(file))
  message(glue("Processing sample: {sample_id}"))
  
  # Lookup cancer type
  cancer_type <- cancer_type_map[[sample_id]]
  if (is.null(cancer_type)) {
    warning(glue("Sample {sample_id} not found in metadata."))
    next
  }
  
  if (!(cancer_type %in% names(data_per_cancer_type))) {
    warning(glue("Cancer type {cancer_type} not found in data_per_cancer_type."))
    next
  }
  
  message(glue("Cancer type: {cancer_type}"))
  
  # Read PHBR CSV
  df <- tryCatch({
    read.csv(file, stringsAsFactors = FALSE)
  }, error = function(e) {
    warning(glue("Could not read {file}: {e$message}"))
    return(NULL)
  })
  
  # Check if valid PHBR file
  if (is.null(df) || nrow(df) == 0) {
    warning(glue("File {basename(file)} is empty or invalid."))
    next
  }
  
  # Add mutation_id column
  df$mutation_id <- paste0(df$gene, "_", gsub("p\\.", "", df$mutation))
  message(glue("Found {length(unique(df$mutation_id))} unique mutations."))
  
  # Get affinity matrix
  aff_mat <- data_per_cancer_type[[cancer_type]]$affinities
  sample_aff_id <- sub("-01$", "", sample_id)
  
  if (!(sample_aff_id %in% rownames(aff_mat))) {
    warning(glue("{sample_aff_id} not in affinity matrix."))
    next
  }
  
  # Map profile scores
  df$prof_score <- sapply(df$mutation_id, function(mut) {
    if (mut %in% colnames(aff_mat)) {
      return(aff_mat[sample_aff_id, mut])
    } else {
      message(glue("Mutation {mut} not found in affinity matrix for {sample_aff_id}"))
      return(NA)
    }
  })
  
  message(glue("Finished scoring {sum(!is.na(df$prof_score))} / {nrow(df)} mutations."))
  
  # Add sample metadata to df
  df$sample_id <- sample_id
  df$cancer_type <- cancer_type
  
  # Store in list
  all_results[[sample_id]] <- df
}

# --------------------------------------------
# Combine results across all samples
# --------------------------------------------
combined_df <- bind_rows(all_results)
message(glue("Combined data from {length(all_results)} samples, total rows: {nrow(combined_df)}"))

# --------------------------------------------
# Filter and prepare data for plotting
# --------------------------------------------
filtered_df <- combined_df %>%
  filter(!is.na(prof_score), !is.na(PHBR_mhc1), !is.na(cancer_type)) %>%
  mutate(
    cancer_type = factor(as.character(cancer_type)),
    label = paste0(gene, ":", gsub("p\\.", "", mutation))
  )

# --------------------------------------------
# Plot PHBR vs profile score
# --------------------------------------------

# Prepare color palette
n_colors <- length(unique(filtered_df$cancer_type))
palette_colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_colors)

# Correlation calculations
pearson <- cor.test(filtered_df$PHBR_mhc1, filtered_df$prof_score, method = "pearson")
spearman <- cor.test(filtered_df$PHBR_mhc1, filtered_df$prof_score, method = "spearman")

# Build plot
ggplot(filtered_df, aes(x = prof_score, y = PHBR_mhc1)) +
  geom_point(aes(color = cancer_type), size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", linetype = "dashed", color = "gray40", se = TRUE, size = 0.6) +
  geom_text_repel(
    aes(label = label),
    size = 2.6,
    box.padding = 0.3,
    max.overlaps = 50,
    segment.color = "gray70",
    segment.size = 0.2
  ) +
  scale_color_manual(values = palette_colors) +
  labs(
    title = "PHBR(MHC-I) Comparison",
    subtitle = sprintf("Pearson r = %.2f, Spearman r = %.2f",
                       pearson$estimate,
                       spearman$estimate),
    x = "PHBR (MHC-I, Claeys et al. 2021)",
    y = "Observed PHBR (MHC-I)",
    color = "Cancer Type",
    caption = "Adapted from Claeys et al., PLOS Genetics 2021 (DOI:10.1371/journal.pgen.1009368)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, margin = margin(b = 10)),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 9, hjust = 0)
  )

# Save plot
ggsave("phbr_vs_profiler.pdf", plot = last_plot(), width = 8, height = 6, units = "in")
