#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ------------------------
# Option parsing
# ------------------------
option_list <- list(
  make_option("--mhc1_files", type = "character", help = "Comma-separated MHC-I .xls files"),
  make_option("--mhc2_files", type = "character", help = "Comma-separated MHC-II .xls files"),
  make_option("--output", type = "character", help = "Output file name"),
  make_option("--kd_threshold", type = "double", default = 500, help = "Threshold for Kd (nM) to define binders [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------
# Helpers
# ------------------------
split_paths <- function(x) {
  gsub("\\[|\\]", "", x) %>%
    strsplit(",") %>%
    .[[1]] %>%
    trimws()
}

mhc1_files <- split_paths(opt$mhc1_files)
mhc2_files <- split_paths(opt$mhc2_files)

mhc1_files <- mhc1_files[file.exists(mhc1_files)]
mhc2_files <- mhc2_files[file.exists(mhc2_files)]

cat("🧪 Valid MHC-I allele files:", length(mhc1_files), "\n")
cat("🧪 Valid MHC-II allele files:", length(mhc2_files), "\n")
cat("📂 Output path:", opt$output, "\n")

if (length(mhc1_files) + length(mhc2_files) == 0) {
  stop("🚫 No valid MHC-I or MHC-II files found.")
}

# ------------------------
# Read one affinity file
# ------------------------
read_affinity_file <- function(file_path) {
  tryCatch({
    cat("📄 Reading:", file_path, "\n")
    df <- read_tsv(file_path, show_col_types = FALSE)

    if (!all(c("nM", "ID", "Rank") %in% colnames(df))) {
      stop("Missing required columns")
    }

    sample_id <- str_extract(basename(file_path), "TCGA-[A-Z0-9]+-[A-Z0-9]+")
    allele <- str_extract(file_path, "HLA-[A-Z0-9:-]+|DRB1_[0-9]+")

    df <- df %>%
      mutate(
        sample_id = sample_id,
        file_id = basename(file_path),  # uniquely tags file
        ID = paste(sample_id, ID, sep = "_"),
        allele = allele
      )

    if ("Peptide" %in% colnames(df)) {
      df <- df %>% mutate(Length = nchar(Peptide))
      if (grepl("mhc1", tolower(file_path))) {
        before <- nrow(df)
        df <- df %>% filter(Length == 9)
        after <- nrow(df)
        cat("🎯 MHC-I: Filtered to 9-mers:", before, "→", after, "\n")
      }
    }

    cat("✅ Loaded:", allele, "| Peptides:", nrow(df), "\n")
    return(df)

  }, error = function(e) {
    message("❌ Failed:", basename(file_path), " — ", e$message)
    return(NULL)
  })
}

# ------------------------
# Robs calculation
# ------------------------
calculate_robs_per_sample <- function(aff_data, class_label = "", kd_threshold = 500) {
  if (is.null(aff_data) || nrow(aff_data) == 0) {
    cat("⚠️ No data for", class_label, "\n")
    return(tibble(
      sample_id = NA_character_,
      n_binding = NA_integer_,
      n_nonbinding = NA_integer_,
      Robs = NA_real_
    ) %>%
      rename_with(~ paste0(., "_", class_label),
                  .cols = c("n_binding", "n_nonbinding", "Robs")))
  }

  harmonic_mean <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0 || any(x == 0)) return(NA_real_)
    return(length(x) / sum(1 / x))
  }

  aff_data %>%
    filter(!is.na(sample_id), !is.na(nM), !is.na(ID), !is.na(file_id)) %>%
    group_by(sample_id, ID, file_id) %>%
    summarise(best_nM = min(nM, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample_id, ID) %>%
    summarise(harmonic_nM = harmonic_mean(best_nM), .groups = "drop") %>%
    mutate(is_binder = harmonic_nM < kd_threshold) %>%
    group_by(sample_id) %>%
    summarise(
      n_binding = sum(is_binder),
      n_nonbinding = sum(!is_binder),
      Robs = ifelse(n_nonbinding == 0, NA, n_binding / n_nonbinding),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(., "_", class_label),
                .cols = c("n_binding", "n_nonbinding", "Robs"))
}

# ------------------------
# Run pipeline
# ------------------------
cat("🔄 Loading MHC-I affinity files...\n")
mhc1_affinities <- map_dfr(mhc1_files, read_affinity_file)

cat("🔄 Loading MHC-II affinity files...\n")
mhc2_affinities <- map_dfr(mhc2_files, read_affinity_file)

cat("📊 Calculating Robs (Kd <", opt$kd_threshold, "nM)...\n")
mhc1_summary <- calculate_robs_per_sample(mhc1_affinities, "MHC1", opt$kd_threshold)
mhc2_summary <- calculate_robs_per_sample(mhc2_affinities, "MHC2", opt$kd_threshold)

final <- full_join(mhc1_summary, mhc2_summary, by = "sample_id")
write_csv(final, opt$output)
cat("✅ Saved Robs summary to:", opt$output, "\n")
