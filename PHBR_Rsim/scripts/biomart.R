#!/usr/bin/env Rscript

# Load the library
library(biomaRt)

# Full output path (corrected to include the directory path)
output_file <- "refseq_to_ensembl.tsv"

# Check if the output file already exists
if (!file.exists(output_file)) {
  cat("🔄 Connecting to Ensembl via biomaRt...\n")
  
  # Connect to Ensembl
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Increase timeout to prevent connection errors
  options(timeout = 300)
  
  # Retrieve the mapping data
  mapping <- getBM(attributes = c("refseq_mrna", "ensembl_transcript_id"), mart = ensembl)
  
  # Clean the data: remove rows with empty or NA refseq_mrna
  mapping <- mapping[!is.na(mapping$refseq_mrna) & mapping$refseq_mrna != "", ]
  
  # Write the cleaned data to the output file
  write.table(mapping, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Confirmation message
  cat("✅ Mapping saved to: ", output_file, "\n")
} else {
  cat("⚠️ Mapping already exists — skipping.\n")
}
