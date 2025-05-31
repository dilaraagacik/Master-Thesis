#!/usr/bin/env Rscript

# === CONFIG ===
args <- commandArgs(trailingOnly = TRUE)
xls_dir <- args[1]
output_file <- args[2]
expected_peptides_per_chunk <- as.integer(ifelse(length(args) >= 3, args[3], 10000))

cat("📁 Scanning directory:", xls_dir, "\n")
cat("📤 Will write:", output_file, "\n")

library(tibble)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# List files
xls_files <- list.files(xls_dir, pattern = "\\.xls$", full.names = TRUE)

# Parse allele + chunk name
parsed <- tibble(
  file = basename(xls_files),
  path = xls_files,
  allele = str_match(basename(xls_files), "(HLA-[A-Z0-9:]+)")[, 1],
  chunk = str_match(basename(xls_files), "part_([a-z]+)\\.xls$")[, 2]
)

# Sort and count
chunk_levels <- sort(unique(parsed$chunk))
parsed <- parsed %>%
  mutate(chunk_index = match(chunk, chunk_levels)) %>%
  arrange(allele, chunk_index)

# Count nM entries
parsed <- parsed %>%
  mutate(valid_nM = map_int(path, function(f) {
    dat <- suppressWarnings(read_tsv(f, skip = 1, col_types = cols(.default = "c")))
    dat %>%
      mutate(nM = suppressWarnings(as.numeric(nM))) %>%
      filter(!is.na(nM)) %>%
      nrow()
  }))

parsed <- parsed %>%
  mutate(start_peptide = (chunk_index - 1) * expected_peptides_per_chunk + 1,
         end_peptide = chunk_index * expected_peptides_per_chunk)

broken_chunks <- parsed %>% filter(valid_nM < expected_peptides_per_chunk)

cat("🔍 Files scanned:", nrow(parsed), "\n")
cat("⚠️ Broken chunks:", nrow(broken_chunks), "\n")

print(broken_chunks %>% select(allele, chunk, valid_nM, start_peptide, end_peptide, file), n = Inf)

write_csv(broken_chunks, output_file)
