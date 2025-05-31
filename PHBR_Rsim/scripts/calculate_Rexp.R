#!/usr/bin/env Rscript

# Suppress unnecessary startup messages and load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(tools)
})

# ---------------------------
# Parse command-line options
# ---------------------------
option_list <- list(
  make_option("--sample_id", type = "character", help = "Sample ID"),
  make_option("--signature_file", type = "character", help = "Signature reference RDS"),
  make_option("--core_rdata", type = "character", help = "Path to GPPM_rand_core.RData"),
  make_option("--mhc1_file", type = "character", help = "MHC-I alleles"),
  make_option("--mhc2_file", type = "character", help = "MHC-II alleles"),
  make_option("--output", type = "character", help = "Output CSV file path"),
  make_option("--metadata_file", type = "character", help = "CSV mapping sample IDs to cancer types")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------
# Detect cancer type for the sample
# ---------------------------
if (is.null(opt$metadata_file)) {
  stop("❌ Metadata file is required to assign cancer type.")
}

# Read metadata and extract cancer type for this sample
metadata <- read.csv(opt$metadata_file, stringsAsFactors = FALSE)
matched <- metadata[metadata$sample_id == opt$sample_id, ]
if (nrow(matched) == 0) stop("❌ Sample ID not found in metadata: ", opt$sample_id)
cancer_type <- matched$cancer_type[1]
cat("📋 Retrieved cancer type from metadata:", cancer_type, "\n")

# ---------------------------
# Load and normalize substitution signature for the cancer type
# ---------------------------
signature_data <- readRDS(opt$signature_file)
signature_cancer <- signature_data %>%
  filter(cancer == cancer_type)

if (nrow(signature_cancer) == 0) {
  stop("❌ No mutation data found for cancer type in signature file: ", cancer_type)
}

# Normalize frequency of trinucleotide substitution types (subst_type3)
ms3 <- signature_cancer %>%
  count(subst_type3) %>%
  filter(!is.na(subst_type3)) %>%
  mutate(freq = n / sum(n)) %>%
  deframe()

if (length(ms3) == 0) {
  stop("❌ No valid subst_type3 values found to create mutation signature for cancer type: ", cancer_type)
}

# ---------------------------
# Read and expand MHC allele lists
# ---------------------------
read_alleles <- function(file_path) {
  lines <- read_lines(file_path)
  alleles <- strsplit(lines, ",") %>% unlist()
  # If allele ends in _twice, include it twice
  unlist(lapply(alleles, function(a) {
    if (grepl("_twice$", a)) rep(gsub("_twice$", "", a), 2) else a
  }))
}

alleles_mhc1 <- read_alleles(opt$mhc1_file)
alleles_mhc2 <- read_alleles(opt$mhc2_file)

if (length(alleles_mhc1) < 2) warning("⚠️ Fewer than 2 MHC-I alleles provided")
if (length(alleles_mhc2) < 2) warning("⚠️ Fewer than 2 MHC-II alleles provided")

# ---------------------------
# Load precomputed simulated peptide data (affinity matrices)
# ---------------------------
cat("📦 Loading simulation data from:", opt$core_rdata, "\n")
sim_vars <- load(opt$core_rdata)
required <- c("subst_type3", "variant", "mhc1_aff", "mhc2_aff")
missing <- setdiff(required, sim_vars)
if (length(missing) > 0) stop("❌ Missing in .RData:", paste(missing, collapse = ", "))

# ---------------------------
# Function to calculate expected immunogenic mutation ratio (Rexp)
# ---------------------------
calculate_Rexp <- function(HLA_genotype, ms3, sim_data, cu = 500) {
  # Calculate harmonic mean for binding prediction across alleles
  harmonic_mean <- function(x) 1 / mean(1 / x, na.rm = TRUE)
  
  # Filter alleles present in simulated affinity matrix
  available_alleles <- HLA_genotype[HLA_genotype %in% colnames(sim_data$HLA_aff_matrix)]
  if (length(available_alleles) < 2) return(NA)
  
  # Compute per-peptide harmonic mean affinity
  gt_matrix <- sim_data$HLA_aff_matrix[, available_alleles, drop = FALSE]
  HLA_aff <- apply(gt_matrix, 1, harmonic_mean)
  isHA <- HLA_aff < cu  # Peptides binding below cutoff
  
  # Build data frame linking substitutions and binding status
  df <- data.frame(
    subst_type3 = sim_data$subst_type3,
    isHA = isHA,
    variant = sim_data$variant
  )
  
  # Create 3D contingency table: substitution x variant type x isHA
  gt_t <- table(df$subst_type3, df$variant, df$isHA)
  
  present_subst <- intersect(rownames(gt_t), names(ms3))
  if (length(present_subst) == 0) return(NA)

  # Apply weights from mutation signature to expected nonsynonymous counts
  sig_vector <- ms3[present_subst]
  weight_vector <- rowSums(gt_t[present_subst, , , drop = FALSE])
  
  safe_weights <- weight_vector > 0 & !is.na(sig_vector)
  if (sum(safe_weights) == 0) return(NA)

  ms_norm <- prop.table(sig_vector[safe_weights] / weight_vector[safe_weights])

  exp_nHLA <- sum(gt_t[present_subst, "nonsynonymous SNV", "TRUE"] * ms_norm, na.rm = TRUE)
  exp_nNonHLA <- sum(gt_t[present_subst, "nonsynonymous SNV", "FALSE"] * ms_norm, na.rm = TRUE)

  if (exp_nNonHLA == 0) return(NA)
  return(exp_nHLA / exp_nNonHLA)  # Expected immunogenicity ratio
}

# ---------------------------
# Compute Rexp for MHC-I
# ---------------------------
sim_data <- list(
  subst_type3 = subst_type3,
  variant = variant,
  HLA_aff_matrix = mhc1_aff
)
Rexp_MHC1 <- calculate_Rexp(alleles_mhc1, ms3, sim_data)

# ---------------------------
# Compute Rexp for MHC-II
# ---------------------------
sim_data$HLA_aff_matrix <- mhc2_aff
Rexp_MHC2 <- calculate_Rexp(alleles_mhc2, ms3, sim_data)

# ---------------------------
# Save results
# ---------------------------
res <- tibble(
  sample_id = opt$sample_id,
  cancer_type = cancer_type,
  Rexp_MHC1 = Rexp_MHC1,
  Rexp_MHC2 = Rexp_MHC2
)

write_csv(res, opt$output)
cat("✅ Rexp values saved to:", opt$output, "\n")
