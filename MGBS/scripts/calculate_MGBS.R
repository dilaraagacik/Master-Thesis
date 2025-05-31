#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(tools)

# Add this just after loading the libraries
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--kd_threshold"), type = "double", default = 500,
              help = "Binding affinity threshold (nM) to define binders [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)

kd_threshold <- opt$options$kd_threshold
args <- opt$args

# Then update the argument order
genotype_csv <- args[1]
mhc1_file <- args[2]
mhc2_file <- args[3]
output_rds <- args[4]

# ----------------------------
# Load input genotype CSV
MHC_gt <- read_csv(genotype_csv, show_col_types = FALSE)

# Ensure unique Patient IDs
if (!"Patient" %in% colnames(MHC_gt)) {
  original_first_col <- colnames(MHC_gt)[1]
  colnames(MHC_gt)[1] <- "Patient"
  message("⚠️ 'Patient' column not found — using first column ('", original_first_col, "') as 'Patient'")
}
MHC_gt$Patient <- make.unique(as.character(MHC_gt$Patient))  # Ensure unique IDs


# Select only relevant genotype columns (A.1, A.2, B.1, B.2, C.1, C.2, DPA1.1, DPA1.2, DPB1.1, DPB1.2)
genotype_columns <- colnames(MHC_gt)
MHC_gt <- MHC_gt[, genotype_columns]
print("genotype_columns")
print(genotype_columns)
# Clean genotype format (convert HLA-A31:01:00 to HLA-A31:01, etc., without asterisk)
clean_hla_genotype <- function(x) {
  # Remove HLA- and convert to proper format (e.g., HLA-A31:01:00 to HLA-A31:01)
  x <- gsub("HLA-", "", x)  # Remove "HLA-" prefix
  x <- sub(":00$", "", x)    # Remove ":00" suffix
  return(x)
}

# Apply the cleaning function to all the genotype columns
for (col in c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2", "DPA1.1", "DPA1.2", "DPB1.1", "DPB1.2")) {
  MHC_gt[[col]] <- clean_hla_genotype(MHC_gt[[col]])
}

# Debug: Print cleaned genotype data
cat("Cleaned Genotype Data:\n")
print(head(MHC_gt))

# ----------------------------
# Load MHC-I and MHC-II affinity matrices (from RDS files)
cat("Loading MHC-I affinity matrix...\n")
MHC1 <- readRDS(mhc1_file)

cat("Loading MHC-II affinity matrix...\n")
MHC2 <- readRDS(mhc2_file)

# Debug: Check the dimensions of MHC1 and MHC2 matrices
cat("Dimensions of MHC1:", dim(MHC1), "\n")
cat("Dimensions of MHC2:", dim(MHC2), "\n")

# Debug: Print first few values of MHC2 (MHC-II affinity matrix)
cat("First few values of MHC2 (MHC-II affinity matrix):\n")
print(head(MHC2))

# ----------------------------
# Define scores (proportion MHC binders) per allele
MGBS1_alleles <- colMeans(MHC1 < kd_threshold)
MGBS2_alleles <- colMeans(MHC2 < kd_threshold)

cat("First few values of MGBS1_alleles:", head(MGBS1_alleles), "\n")
#cat("First few values of MGBS2_alleles:", head(MGBS2_alleles), "\n")

# ----------------------------
# MHC-I genotypes in the same format as affinities
for (a in c("A", "B", "C")) {
  for (i in 1:2) {
    col <- paste0(a, ".", i)
    MHC_gt[[col]] <- ifelse(
      is.na(MHC_gt[[col]]),
      NA,
      paste0("HLA-", a, MHC_gt[[col]])
    )
  }
}


MHC1_gt <- MHC_gt[, 2:7]  # Only the HLA allele columns for MHC-I

cat("First few rows of MHC1 genotypes (MHC1_gt):\n")
print(head(MHC1_gt))

# Debug: Check which alleles are missing from affinity matrix
missing_alleles <- setdiff(unique(unlist(MHC1_gt)), names(MGBS1_alleles))
cat("⚠️ Alleles in genotype data but missing in affinity matrix:\n")
print(missing_alleles)

# Count how many alleles are missing per patient
MGBS1_na_counts <- apply(MHC1_gt, 1, function(x) sum(is.na(MGBS1_alleles[as.character(x)])))

# Get indices to keep: less than 3 missing alleles
valid_indices <- which(MGBS1_na_counts < 3)

# Filter genotype and allele data
MHC_gt <- MHC_gt[valid_indices, ]
MHC1_gt <- MHC1_gt[valid_indices, , drop = FALSE]

cat("✅ Removed", length(MGBS1_na_counts) - length(valid_indices), "patients with ≥3 missing alleles\n")


# Calculate MGBS-I (Mean of MGBS1 alleles for each genotype)
MGBS1 <- apply(MHC1_gt, 1, function(x) mean(MGBS1_alleles[as.character(x)], na.rm = TRUE))
names(MGBS1) <- MHC_gt$Patient

cat("First few values of MGBS1:", head(MGBS1), "\n")

# ----------------------------
# MHC-II genotypes in the same format as affinities (heterodimers)
# Clean MHC-II genotypes (remove colons and trailing ":00")
for (col in grep("D", colnames(MHC_gt), value = TRUE)) {
  MHC_gt[[col]] <- gsub(":", "", MHC_gt[[col]])  # Remove all colons
  MHC_gt[[col]] <- sub("00$", "", MHC_gt[[col]])  # Remove trailing 00
}
MHC2_gt <- data.frame(
  row.names = MHC_gt$Patient,
  HLA_DP_11 = paste0("HLA-", "DPA1", MHC_gt[["DPA1.1"]], "-", "DPB1", MHC_gt[["DPB1.1"]]),
  HLA_DP_12 = paste0("HLA-", "DPA1", MHC_gt[["DPA1.1"]], "-", "DPB1", MHC_gt[["DPB1.2"]]),
  HLA_DP_21 = paste0("HLA-", "DPA1", MHC_gt[["DPA1.2"]], "-", "DPB1", MHC_gt[["DPB1.1"]]),
  HLA_DP_22 = paste0("HLA-", "DPA1", MHC_gt[["DPA1.2"]], "-", "DPB1", MHC_gt[["DPB1.2"]]),
  HLA_DQ_11 = paste0("HLA-", "DQA1", MHC_gt[["DQA1.1"]], "-", "DQB1", MHC_gt[["DQB1.1"]]),
  HLA_DQ_12 = paste0("HLA-", "DQA1", MHC_gt[["DQA1.1"]], "-", "DQB1", MHC_gt[["DQB1.2"]]),
  HLA_DQ_21 = paste0("HLA-", "DQA1", MHC_gt[["DQA1.2"]], "-", "DQB1", MHC_gt[["DQB1.1"]]),
  HLA_DQ_22 = paste0("HLA-", "DQA1", MHC_gt[["DQA1.2"]], "-", "DQB1", MHC_gt[["DQB1.2"]]),
  HLA_DR_1 = paste0("DRB1_", MHC_gt[["DRB1.1"]]),
  HLA_DR_2 = paste0("DRB1_", MHC_gt[["DRB1.2"]])
)

# Calculate MGBS-II (Mean of MGBS2 alleles for each genotype)
MGBS2 <- apply(MHC2_gt, 1, function(x) mean(MGBS2_alleles[as.character(x)], na.rm = TRUE))

cat("First few values of MGBS2:", head(MGBS2), "\n")
cat("First few rows of MHC2 genotypes (MHC2_gt):\n")
print(head(MHC2_gt))
# ----------------------------
# Check for NA values in MGBS1 and MGBS2
cat("NA values in MGBS1:", sum(is.na(MGBS1)), "\n")
cat("NA values in MGBS2:", sum(is.na(MGBS2)), "\n")

# ----------------------------
# Ensure that MGBS2 is aligned with MGBS1 (based on patient names)
MGBS2 <- MGBS2[names(MGBS1)]  # Align MGBS2 by the names (Patient IDs) of MGBS1

cat("After alignment, First few values of MGBS2:", head(MGBS2), "\n")
# Save MGBS1 and MGBS2 separately
MGBS1_df <- data.frame(
  patient = as.character(names(MGBS1)),
  MGBS1 = as.numeric(MGBS1),
  stringsAsFactors = FALSE
)

MGBS2_df <- data.frame(
  patient = as.character(names(MGBS2)),
  MGBS2 = as.numeric(MGBS2),
  stringsAsFactors = FALSE
)

saveRDS(MGBS1_df, "MGBS1_scores.rds")
saveRDS(MGBS2_df, "MGBS2_scores.rds")
