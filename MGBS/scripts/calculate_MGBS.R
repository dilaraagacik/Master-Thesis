#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(tools)
suppressPackageStartupMessages(library(optparse))

# Parse options
option_list <- list(
  make_option(c("--kd_threshold"), type = "double", default = 500,
              help = "Binding affinity threshold (nM) to define binders [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
kd_threshold <- opt$options$kd_threshold
args <- opt$args

# Argument order
genotype_csv <- args[1]
mhc1_file <- args[2]
mhc2_file <- args[3]
output_rds <- args[4]

# ----------------------------
# Load input genotype CSV
MHC_gt <- read_csv(genotype_csv, col_types = cols(.default = "c"))

# Ensure unique Patient IDs
if (!"Patient" %in% colnames(MHC_gt)) {
  original_first_col <- colnames(MHC_gt)[1]
  colnames(MHC_gt)[1] <- "Patient"
  message("'Patient' column not found — using first column ('", original_first_col, "') as 'Patient'")
}
MHC_gt$Patient <- make.unique(as.character(MHC_gt$Patient))

# Debug
cat("DEBUG — column names in MHC_gt:\n")
print(colnames(MHC_gt))
cat("DEBUG — first few rows of MHC_gt:\n")
print(head(MHC_gt))
cat("DEBUG — parsing problems:\n")
print(problems(MHC_gt))

# Clean genotype format
clean_hla_genotype <- function(x) {
  x <- gsub("HLA-", "", x)
  x <- sub(":00$", "", x)
  return(x)
}
for (col in c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2", "DPA1.1", "DPA1.2", "DPB1.1", "DPB1.2")) {
  MHC_gt[[col]] <- clean_hla_genotype(MHC_gt[[col]])
}

# ----------------------------
# Load MHC-I and MHC-II affinity matrices
cat("Loading MHC-I affinity matrix...\n")
MHC1 <- readRDS(mhc1_file)

cat("Loading MHC-II affinity matrix...\n")
MHC2 <- readRDS(mhc2_file)

# ----------------------------
# Define scores (proportion MHC binders) per allele
MGBS1_alleles <- colMeans(MHC1 < kd_threshold)
MGBS2_alleles <- colMeans(MHC2 < kd_threshold)

# ----------------------------
# MHC-I genotypes in same format as affinities
for (a in c("A", "B", "C")) {
  for (i in 1:2) {
    col <- paste0(a, ".", i)
    MHC_gt[[col]] <- paste0("HLA-", a, MHC_gt[[col]])
  }
}
rownames(MHC_gt) <- MHC_gt$Patient
MHC1_gt <- MHC_gt[, 2:7]

# Check missing alleles — MHC-I
missing_alleles_MHC1 <- setdiff(unique(unlist(MHC1_gt)), names(MGBS1_alleles))
cat("Missing alleles in MHC-I affinity matrix:\n")
print(missing_alleles_MHC1)

# Filter patients based on missing alleles (Claeys logic)
MGBS1_na_counts <- apply(MHC1_gt, 1, function(x) sum(is.na(MGBS1_alleles[as.character(x)])))
valid_indices1 <- which(MGBS1_na_counts < 3)
cat("Keeping", length(valid_indices1), "patients for MHC-I (removed", nrow(MHC_gt) - length(valid_indices1), ")\n")

# ----------------------------
# MHC-II genotypes in same format as affinities
for (col in grep("D", colnames(MHC_gt), value = TRUE)) {
  MHC_gt[[col]] <- gsub(":", "", MHC_gt[[col]])
  MHC_gt[[col]] <- sub("00$", "", MHC_gt[[col]])
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

# Check missing alleles — MHC-II
missing_alleles_MHC2 <- setdiff(unique(unlist(MHC2_gt)), names(MGBS2_alleles))
cat("Missing alleles in MHC-II affinity matrix:\n")
print(missing_alleles_MHC2)

# Filter patients based on missing alleles — MHC-II
MGBS2_na_counts <- apply(MHC2_gt, 1, function(x) sum(is.na(MGBS2_alleles[as.character(x)])))
valid_indices2 <- which(MGBS2_na_counts < 3)
cat("Keeping", length(valid_indices2), "patients for MHC-II (removed", nrow(MHC_gt) - length(valid_indices2), ")\n")

# ----------------------------
# Final intersection — patients valid for BOTH MHC-I and MHC-II
final_indices <- intersect(valid_indices1, valid_indices2)
cat("Final patients kept:", length(final_indices), "/", nrow(MHC_gt), "\n")

# Filter data
MHC_gt <- MHC_gt[final_indices, ]
MHC1_gt <- MHC1_gt[final_indices, , drop = FALSE]
MHC2_gt <- MHC2_gt[final_indices, , drop = FALSE]

rownames(MHC1_gt) <- MHC_gt$Patient
rownames(MHC2_gt) <- MHC_gt$Patient

# ----------------------------
# Calculate MGBS1 and MGBS2 for final patients
MGBS1 <- apply(MHC1_gt, 1, function(x) mean(MGBS1_alleles[as.character(x)], na.rm = TRUE))
MGBS2 <- apply(MHC2_gt, 1, function(x) mean(MGBS2_alleles[as.character(x)], na.rm = TRUE))

# Set names explicitly for safety
names(MGBS1) <- rownames(MHC1_gt)
names(MGBS2) <- rownames(MHC2_gt)

# ----------------------------
# Merge results
MGBS_df <- data.frame(
  patient = names(MGBS1),
  MGBS1 = MGBS1,
  MGBS2 = MGBS2[names(MGBS1)]
)

keep_idx <- which(!is.na(MGBS_df$MGBS1) & !is.na(MGBS_df$MGBS2))
MGBS_df <- MGBS_df[keep_idx, ]
cat("Final patients kept for ECDF calculation:", nrow(MGBS_df), "\n")

# Add MGBSd
calculate_MGBS_norm <- function(df, var_MHC1, var_MHC2) {
  m1 <- df[, var_MHC1]
  m2 <- df[, var_MHC2]
  
  m1f <- ecdf(discard(m1, is.na))
  m2f <- ecdf(discard(m2, is.na))
  
  mdnorm <- set_names(m1f(m1) - m2f(m2), rownames(df))
  
  return(mdnorm)
}
MGBS_df$MGBSd <- calculate_MGBS_norm(MGBS_df, "MGBS1", "MGBS2")

# ----------------------------
# Save full dataframe
cat("Saving MGBS_df to:", output_rds, "\n")
saveRDS(MGBS_df, output_rds)
