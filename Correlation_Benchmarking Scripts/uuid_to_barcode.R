library(data.table)
library(tidyverse)

# Init
star_dir <- ".../downloads/TCGA/mRNA/STAR_counts"
uuid_list <- list.files(star_dir)

# For testing → small subset first!
uuid_list_subset <- uuid_list[1:10]  # Take first 10 for test

#First: read one file to get the gene list
first_uuid <- uuid_list_subset[1]
first_file <- list.files(file.path(star_dir, first_uuid), pattern = "\\.tsv$", full.names = TRUE)[1]

df_first <- fread(first_file)

# Use gene_name instead of gene_id
gene_names <- df_first$gene_name
num_genes <- length(gene_names)

# Initialize big matrix
big_mat <- matrix(NA, nrow = num_genes, ncol = length(uuid_list_subset))
rownames(big_mat) <- gene_names
colnames(big_mat) <- uuid_list_subset

#  Loop with progress bar
pb <- txtProgressBar(min = 0, max = length(uuid_list_subset), style = 3)

for (i in seq_along(uuid_list_subset)) {
  uuid <- uuid_list_subset[i]
  file_path <- list.files(file.path(star_dir, uuid), pattern = "\\.tsv$|\\.txt$", full.names = TRUE)[1]
  
  if (is.na(file_path) || !file.exists(file_path)) {
    message(" Missing file for UUID: ", uuid)
    next
  }
  
  df <- fread(file_path)
  
  # No grep — keep all rows, match on gene_name
  big_mat[, i] <- df$unstranded[match(gene_names, df$gene_name)]
  
  setTxtProgressBar(pb, i)
}
close(pb)

#  Save test matrix
saveRDS(big_mat, file = "tcga_gene_counts_matrix_gene_name_subset_test.rds")
message("Test matrix saved: tcga_gene_counts_matrix_gene_name_subset_test.rds")

