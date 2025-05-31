#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript parse_netmhc_to_rds.R <xls_dir> <output_rds_file>")
}

xls_dir <- args[1]
output_file <- args[2]

files <- list.files(path = xls_dir, pattern = "\\.xls$", full.names = TRUE)
if (length(files) == 0 || all(grepl("dummy", files))) {
  cat("⚠️ No real NetMHC output files found. Skipping RDS creation.\n")
  quit(save = "no", status = 0)
}

cat("📦 Found", length(files), "NetMHC output files\n")

# Function to extract allele name from filename (e.g. HLA-A11:01)
extract_allele <- function(file) {
  match <- str_match(basename(file), "HLA-\\w+\\d+:\\d+")
  if (is.na(match[1])) return("unknown")
  return(match[1])
}

# Group files by allele
allele_names <- map_chr(files, extract_allele)
cat("🧬 Unique alleles extracted:", paste(unique(allele_names), collapse = ", "), "\n")
allele_groups <- split(files, allele_names)

# For each allele, read all its chunks and combine
affinity_list <- imap(allele_groups, function(file_group, allele) {
  cat("\n📥 Processing allele:", allele, "with", length(file_group), "files\n")
  
  dfs <- map(file_group, function(f) {
    cat("  📂 Reading file:", f, "\n")
    dat <- read_tsv(f, skip = 1, col_types = cols(.default = "c"))
    cat("  ✅ Read", nrow(dat), "rows\n")
    
    if (!"nM" %in% colnames(dat)) {
      cat("  ⚠️ Column 'nM' not found in:", f, "\n")
      return(tibble(nM = numeric(0)))
    }
    
    dat <- dat %>% filter(!is.na(nM)) %>% transmute(nM = as.numeric(nM))
    cat("  📊 Filtered to", nrow(dat), "valid 'nM' rows\n")
    return(dat)
  })
  
  all_nM <- bind_rows(dfs)
  cat("📊 Combined rows for allele", allele, ":", nrow(all_nM), "\n")
  tibble(!!sym(allele) := all_nM$nM)
})

# Build full matrix: rows = peptides (assumed to be same order), cols = alleles
# Pad columns to the same number of rows
max_len <- max(map_int(affinity_list, nrow))
affinity_list <- map(affinity_list, function(df) {
  if (nrow(df) < max_len) {
    df[(nrow(df)+1):max_len, ] <- NA  # pad with NAs
  }
  df
})

aff_matrix <- bind_cols(affinity_list)


saveRDS(aff_matrix, output_file)
cat("✅ Saved affinity matrix to:", output_file, "\n")
