#!/usr/bin/env Rscript

library(jsonlite)
library(dplyr)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
mapping_file <- args[1]
output_file <- args[2]
json_files <- args[3:length(args)]

# Filter out missing or bad paths
json_files <- json_files[!is.na(json_files) & file.exists(json_files)]

cat("📦 Using the following JSON files:\n")
print(json_files)

if (length(json_files) == 0) {
  stop("❌ No valid JSON files were passed to the script.")
}

# Load and prepare mapping
mapping <- read.csv(mapping_file, stringsAsFactors = FALSE)
mapping$tissue <- ifelse(grepl("metastasized", mapping$tissue, ignore.case = TRUE), "M",
                         ifelse(grepl("Primary", mapping$tissue, ignore.case = TRUE), "P", "N"))
mapping$patient <- paste(mapping$subject, mapping$tissue, sep = "_")
mapping$Run <- trimws(mapping$Run)

# Function to extract the XX:YY part of HLA alleles
parse_allele <- function(alleles) {
  str_extract(alleles, "\\d{2}:\\d{2}")
}

# Parse genotype files
raw_data <- list()

for (file in json_files) {
  genotype <- fromJSON(file)
  
  # Extract SRR from filename
  srr <- str_extract(basename(file), "SRR[0-9]+")
  if (length(srr) != 1 || is.na(srr)) next

  row <- c(
    A.1     = parse_allele(genotype$A)[1],
    A.2     = parse_allele(genotype$A)[2],
    B.1     = parse_allele(genotype$B)[1],
    B.2     = parse_allele(genotype$B)[2],
    C.1     = parse_allele(genotype$C)[1],
    C.2     = parse_allele(genotype$C)[2],
    DPA1.1  = parse_allele(genotype$DPA1)[1],
    DPA1.2  = parse_allele(genotype$DPA1)[2],
    DPB1.1  = parse_allele(genotype$DPB1)[1],
    DPB1.2  = parse_allele(genotype$DPB1)[2],
    DRB1.1  = parse_allele(genotype$DRB1)[1],  # Add DRB1
    DRB1.2  = parse_allele(genotype$DRB1)[2],  # Add DRB1
    DQA1.1  = parse_allele(genotype$DQA1)[1],  # Add DQA1
    DQA1.2  = parse_allele(genotype$DQA1)[2],  # Add DQA1
    DQB1.1  = parse_allele(genotype$DQB1)[1],  # Add DQB1
    DQB1.2  = parse_allele(genotype$DQB1)[2]   # Add DQB1
  )

  raw_data[[srr]] <- row
}


allele_df <- as.data.frame(do.call(rbind, raw_data), stringsAsFactors = FALSE)
allele_df <- tibble::rownames_to_column(allele_df, "SRR")

# Ensure unique row names for SRR IDs
if (any(duplicated(allele_df$SRR))) {
  allele_df$SRR <- make.unique(allele_df$SRR)  # Adds suffix to duplicates
}

# Assign row names
rownames(allele_df) <- allele_df$SRR


cat("🔎 Parsed allele dataframe:\n")
print(allele_df)

# Merge with mapping
final_df <- merge(mapping[, c("Run", "patient")], allele_df, by.x = "Run", by.y = "SRR")

if (nrow(final_df) == 0) {
  stop("❌ No matching SRR IDs found in SraRunTable.csv for the parsed genotype files.")
}

final_df <- final_df[, -1, drop = FALSE]
colnames(final_df)[1] <- "Patient"

cat("✅ Final patient-allele dataframe:\n")
print(final_df)

# Save to RDS
write.csv(final_df, file = paste0(tools::file_path_sans_ext(output_file), ".csv"), row.names = FALSE)
cat("✅ Saved output to:", output_file, "\n")
