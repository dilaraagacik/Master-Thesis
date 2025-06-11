library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ConsensusTME)
library(edgeR)
library(scales)
library(pheatmap)
library(broom)
library(tidyr)


# Read file (first column = gene, other columns = samples)
expr <- readRDS("~/Desktop/tcga_gene_counts_matrix_with_barcodes.rds")
meta_data = readRDS("~/Desktop/patient_info_from_uuid_SERVER_ALL_tcga.rds")
barcode_to_cancer <- setNames(meta_data$cancer_type, meta_data$submitter_id)
colnames(expr) <- sub("\\.\\.\\..*", "", colnames(expr))

length(unique(na.omit(colnames(expr))))
length(unique(na.omit(meta_data$submitter_id)))


# Check the structure: usually, first column is Geneid/Ensembl and the rest are sample columns
head(expr)[, 1:6]
ncol(expr) # Check number of columns
# Assume expr is your tibble
expr_df <- as.data.frame(expr)
rownames(expr_df) <- expr_df[[1]]
clean_names <- function(x) {
  gsub("\\.[0-9]+$", "", x)
}

# 1. Replace NA with 0
expr_df_clean <- expr_df
expr_df_clean[is.na(expr_df_clean)] <- 0

# 2. Remove duplicated columns, keep first
dup_cols <- duplicated(colnames(expr_df_clean))
expr_df_clean <- expr_df_clean[, !dup_cols]

# 3. Remove columns with library size = 0
lib_sizes <- colSums(expr_df_clean)
expr_df_clean <- expr_df_clean[, lib_sizes > 0]

# 4. Now safe to create DGEList
dge <- DGEList(counts = expr_df_clean)

# 5. Now safe to compute CPM
cpm_matrix <- cpm(dge, log = TRUE)




# Make sure sample_cancer is a named vector: names are sample IDs, values are cancer types
# E.g., sample_cancer <- setNames(meta$cancer_type, meta$barcode)
unique_cancers <- unique(na.omit(meta_data$cancer_type))
sample_cancer  <- setNames(meta_data$cancer_type, meta_data$submitter_id)


cancer_results <- list()
batch_size <- 200


for (cancer in unique_cancers) {
  cat("Processing cancer type:", cancer, "\n")
  cancer_samples <- names(sample_cancer[sample_cancer == cancer])
  cancer_samples <- intersect(cancer_samples, colnames(cpm_matrix))
  
  if (length(cancer_samples) == 0) {
    warning("No matching samples found for cancer type: ", cancer)
    next
  }
  
  cancer_sample_batches <- split(cancer_samples, ceiling(seq_along(cancer_samples) / batch_size))
  cancer_batch_results <- list()
  
  for (i in seq_along(cancer_sample_batches)) {
    cat("  Batch", i, "of", length(cancer_sample_batches), "for", cancer, "\n")
    batch_samples <- cancer_sample_batches[[i]]
    batch_mat <- as.matrix(cpm_matrix[, batch_samples, drop = FALSE])
    mode(batch_mat) <- "numeric"
    
    result <- tryCatch({
      ConsensusTME::consensusTMEAnalysis(
        batch_mat, 
        cancer = if (cancer %in% ConsensusTME::cancerAll) cancer else "unfiltered", 
        statMethod = "ssgsea"
      )
    }, error = function(e) {
      warning("    ConsensusTME failed for cancer type: ", cancer, " (batch ", i, "). Skipping batch.")
      NULL
    })
    
    if (!is.null(result)) {
      cancer_batch_results[[i]] <- result
    }
  }
  
  if (length(cancer_batch_results) > 0) {
    cancer_results[[cancer]] <- do.call(cbind, cancer_batch_results)
    cat("  --> Processed", ncol(cancer_results[[cancer]]), "samples for", cancer, "\n")
  } else {
    warning("No successful batches for cancer type: ", cancer)
  }
}

#Total patient who have immune deconvolution results
total_patients <- sum(sapply(cancer_results, ncol))
print(total_patients)

#

# Create data frame from named vector
df <- data.frame(
  Patient = names(sample_cancer),
  Cancer = sample_cancer,
  stringsAsFactors = FALSE
)

# Remove duplicate entries per patient per cancer type
df_unique <- df[!duplicated(df), ]

# Count unique patients per cancer type
unique_counts <- table(df_unique$Cancer)

# View result
print(unique_counts)
sum(unique_counts)


# Get all immune signatures
immune_sigs <- rownames(cancer_results[[1]])

# Initialize matrix: rows = immune cell types, columns = cancer types
summary_mat <- matrix(NA, nrow = length(immune_sigs), ncol = length(cancer_results))
rownames(summary_mat) <- immune_sigs
colnames(summary_mat) <- names(cancer_results)

# Fill in the matrix
for (cancer in names(cancer_results)) {
  mat <- cancer_results[[cancer]]
  common_sigs <- intersect(rownames(mat), immune_sigs)
  summary_mat[common_sigs, cancer] <- apply(mat[common_sigs, , drop=FALSE], 1, median, na.rm=TRUE)
}


summary_mat_scaled <- t(scale(t(summary_mat)))
rownames(summary_mat_scaled ) <- gsub("_", " ", rownames(summary_mat_scaled ))

set.seed(42) 
pheatmap(
  summary_mat_scaled,
  colors <- colorRampPalette(c("#053061", "#2166ac", "#4393c3", "#92c5de", 
                               "#d1e5f0", "#fddbc7", "#f4a582", "#d6604d", 
                               "#b2182b", "#67001f"))(100)
  ,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  filename = "ImmuneDeconv_CancerSummary.png",
  main ="Immune Cell Landscape Across Different Cancer Types in the TCGA Dataset",
  width = 10, height = 8
)

# Save to file
saveRDS(cancer_results , "ConsensusTME_results_geneCounts.rds")
write.csv(final_result, "ConsensusTME_results.csv")

library(tibble)

summary_mat_scaled[ , "MESO"]

# Convert to data frame
meso_summary <- as.data.frame(summary_mat_scaled[ , "MESO"]) %>%
  rownames_to_column("Immune_Cell") %>%
  rename(MESO_Score = `summary_mat_scaled[, "MESO"]`)

ggplot(meso_summary, aes(x = reorder(Immune_Cell, MESO_Score), y = MESO_Score, fill = MESO_Score)) +
  geom_col(width = 0.6, color = "grey30", alpha = 0.95) +
  coord_flip() +
  scale_fill_gradient(low = "#92c5de", high = "#b2182b") +
  theme_minimal(base_size = 16) +
  labs(
    title = "Immune Cell Infiltration Summary — MESO",
    x = "",
    y = "Scaled Immune Score",
    fill = "Score"
  ) +
  theme(
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    legend.position = "right"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
ggsave("ImmuneCell_Infiltration_MESO.png", width = 10, height = 6)
