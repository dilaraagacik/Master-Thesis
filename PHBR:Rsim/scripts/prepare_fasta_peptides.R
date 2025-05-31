#!/usr/bin/env Rscript

library(dplyr)
library(readr)

# --- Handle input arguments ---
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1:(length(args)-2)]
mhc1_out <- args[length(args)-1]
mhc2_out <- args[length(args)]

# --- Load peptide data ---
message("🔍 Reading input peptide tables...")
peptides <- bind_rows(lapply(input_files, read.delim, stringsAsFactors = FALSE))

# --- Prepare FASTA headers ---
peptides <- peptides %>%
  mutate(header = paste0(">", Gene, "_", AAChange)) %>%
  filter(!is.na(Peptide_MHC1) | !is.na(Peptide_MHC2)) %>%
  distinct(header, Peptide_MHC1, Peptide_MHC2, .keep_all = TRUE)

# --- Write MHC-I FASTA ---
message("📤 Writing MHC-I FASTA: ", mhc1_out)
mhc1_fasta <- unlist(apply(peptides, 1, function(row) {
  c(row["header"], row["Peptide_MHC1"])
}))
writeLines(mhc1_fasta, mhc1_out)

# --- Write MHC-II FASTA ---
message("📤 Writing MHC-II FASTA: ", mhc2_out)
mhc2_fasta <- unlist(apply(peptides, 1, function(row) {
  c(row["header"], row["Peptide_MHC2"])
}))
writeLines(mhc2_fasta, mhc2_out)

message("✅ FASTA preparation complete.")
