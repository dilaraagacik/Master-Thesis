#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)

# Arguments
counts_file <- args[1]
output_folder <- args[2]
mapping_file <- args[3]

message("🔄 Reading input files...")
counts <- read_tsv(counts_file)
genes_mapping <- read_csv(mapping_file)


message("🔍 Merging gene mapping...")
genes_mapping_filtered <- genes_mapping %>%
  filter(ensembl_gene_id %in% counts$Geneid)


counts$ensembl_gene_id <- counts$Geneid
counts_mapped <- inner_join(counts,genes_mapping_filtered, by = "ensembl_gene_id")


message("📉 Initial gene count: ", nrow(counts))
message("📈 Genes after mapping: ", nrow(counts_mapped))

# Remove duplicates
counts_mapped <- counts_mapped %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "") %>%  # remove blank gene names
  distinct(hgnc_symbol, .keep_all = TRUE)




counts_ready <- counts_mapped %>%
  column_to_rownames("hgnc_symbol") %>%
  select(-Geneid, -ensembl_gene_id)

message("🧪 Checking for zero-expression genes per patient...")

# Create output folder
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Write sample metadata
message("🧾 Sample IDs in counts: ", paste(head(colnames(counts_ready)), collapse = ", "), " ...")

# Save a debug snapshot of first few rows
write_csv(head(counts_ready), file.path(output_folder, "DEBUG_counts_ready_head.csv"))

# For each patient/sample, find zero-expression genes
for (patient in colnames(counts_ready)) {
  expr_values <- counts_ready[[patient]]
  zero_genes <- rownames(counts_ready)[expr_values == 0]
  
  message("🧬 Sample: ", patient, " → Zero genes: ", length(zero_genes))
  
  output_path <- file.path(output_folder, paste0(patient, "_zero_genes.txt"))
  writeLines(zero_genes, output_path)
}
