#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
  library(stringr)
})
# Usage: Rscript hla_parser.R output_directory > filtered_hla.csv

args <- commandArgs(trailingOnly = TRUE)
input_dir <- ifelse(length(args) == 0, ".", args[1])

# List all *_HLAgenotype.txt files
files <- list.files(input_dir, pattern = "_HLAgenotype\\.txt$", full.names = TRUE)

# Define HLA genes of interest
genes <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")

# Prepare header
header <- c("sample_id", unlist(lapply(genes, function(g) c(paste0(g, ".1"), paste0(g, ".2")))))
cat(paste(header, collapse = ","), "\n")

# Process each file
for (file in files) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_id <- sub("_HLAgenotype\\.txt$", "", basename(file))
  
  values <- c(sample_id)
  for (g in genes) {
    row <- df[df$Gene == g, ]
    if (nrow(row) == 1) {
      values <- c(values, row$Allele1, row$Allele2)
    } else {
      values <- c(values, "NA", "NA")
    }
  }
  cat(paste(values, collapse = ","), "\n")
}
