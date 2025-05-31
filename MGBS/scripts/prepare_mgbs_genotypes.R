#!/usr/bin/env Rscript

library(jsonlite)
library(dplyr)
library(stringr)
library(tibble)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
output_file <- args[1]
json_files <- args[-1]

# Validate input files
json_files <- json_files[!is.na(json_files) & file.exists(json_files)]
if (length(json_files) == 0) stop("❌ No valid JSON files provided.")

cat("📦 Processing the following JSON files:\n")
print(json_files)

# Helper functions
parse_allele <- function(alleles) str_extract(alleles, "\\d{2}:\\d{2}")
extract_sample_id <- function(file) str_replace(basename(file), "_filtered_genotype\\.json$", "")

# Read and parse genotype JSONs
allele_data <- lapply(json_files, function(file) {
  g <- fromJSON(file)
  sample_id <- extract_sample_id(file)

  data.frame(
    Sample  = sample_id,
    A.1     = parse_allele(g$A)[1],
    A.2     = parse_allele(g$A)[2],
    B.1     = parse_allele(g$B)[1],
    B.2     = parse_allele(g$B)[2],
    C.1     = parse_allele(g$C)[1],
    C.2     = parse_allele(g$C)[2],
    DPA1.1  = parse_allele(g$DPA1)[1],
    DPA1.2  = parse_allele(g$DPA1)[2],
    DPB1.1  = parse_allele(g$DPB1)[1],
    DPB1.2  = parse_allele(g$DPB1)[2],
    DQA1.1  = parse_allele(g$DQA1)[1],
    DQA1.2  = parse_allele(g$DQA1)[2],
    DQB1.1  = parse_allele(g$DQB1)[1],
    DQB1.2  = parse_allele(g$DQB1)[2],
    DRB1.1  = parse_allele(g$DRB1)[1],
    DRB1.2  = parse_allele(g$DRB1)[2],
    stringsAsFactors = FALSE
  )
})

final_df <- bind_rows(allele_data)

cat("✅ Final sample-allele genotype table:\n")
print(final_df)

# Save to CSV
write.csv(final_df, file = paste0(file_path_sans_ext(output_file), ".csv"), row.names = FALSE)
cat("✅ Saved output to:", output_file, "\n")
