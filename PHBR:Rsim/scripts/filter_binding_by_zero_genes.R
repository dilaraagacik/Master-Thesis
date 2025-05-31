#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
  library(stringr)
})

option_list = list(
  make_option("--mhc1_files", type = "character"),
  make_option("--mhc2_files", type = "character"),
  make_option("--zero_genes", type = "character"),
  make_option("--out_prefix", type = "character")
)

opt = parse_args(OptionParser(option_list = option_list))

mhc1_files = strsplit(opt$mhc1_files, ",")[[1]]
mhc2_files = strsplit(opt$mhc2_files, ",")[[1]]
zero_genes_file = opt$zero_genes
out_prefix = opt$out_prefix

genes_to_exclude = readLines(zero_genes_file)

filter_binding_file <- function(files, type) {
  for (file in files) {
    cat("\n📂 Processing file:", file, "\n")

    if (!file.exists(file)) {
      warning(paste("⚠️ File missing:", file))
      next
    }

    # Load NetMHC output
    df <- tryCatch(read.delim(file, row.names = NULL, header = TRUE, skip = 1), error = function(e) {
      warning(paste("❌ Failed to read:", file, "->", e$message))
      return(NULL)
    })

    if (is.null(df) || !"ID" %in% colnames(df)) {
      warning(paste("⚠️ Skipping", file, "- no 'ID' column."))
      next
    }

    if (!"Peptide" %in% colnames(df)) {
      warning("⚠️ Missing 'Peptide' column — cannot filter by length.")
      df$Length <- NA
    } else {
      df$Length <- nchar(df$Peptide)
    }

    # Extract gene from ID
    df$gene <- str_extract(df$ID, "^[^_]+")

    # Extract sample ID from filename
    sample_id <- str_extract(file, "TCGA-[A-Z0-9-]+")
    if (is.na(sample_id)) {
      warning(paste("⚠️ Could not find sample ID in", file))
      sample_id <- "UnknownSample"
    }

    df$sample_id <- sample_id
    zero_gene_file <- file.path(dirname(opt$zero_genes), paste0(sample_id, "_zero_genes.txt"))

    if (!file.exists(zero_gene_file)) {
      warning("⚠️ No zero-expression gene file found for", sample_id, "→ skipping expression filter.")
      genes_to_exclude <- character(0)
    } else {
      genes_to_exclude <- readLines(zero_gene_file)
      cat("🧬 Zero-expression genes loaded for", sample_id, ":", length(genes_to_exclude), "\n")
    }

    # Track pre-filter stats
    before_n <- nrow(df)
    before_zero_matches <- sum(df$gene %in% genes_to_exclude)

    # Apply filtering
    df_filtered <- df %>%
      filter(!(gene %in% genes_to_exclude)) %>%
      filter((type == "filtered_mhc1" & Length == 9) | (type == "filtered_mhc2" & Length == 15))

    after_n <- nrow(df_filtered)
    cat(sprintf("✅ Filtered: %s | Before: %d → After: %d | Removed due to zero-expression: %d\n",
                basename(file), before_n, after_n, before_zero_matches))

    if (after_n == 0) {
      warning("⚠️ All peptides removed after filtering — check gene symbols and input files.")
    }

    allele <- str_extract(basename(file), "HLA-[A-Z0-9:-]+|DRB1_[0-9]+")
    out_file <- paste0(type, "_", sample_id, "_", allele, ifelse(grepl("_twice", file), "_twice", ""), ".xls")

    write_delim(df_filtered, out_file, delim = "\t")
    cat("💾 Saved:", out_file, "\n")
  }
}


# Apply filtering
filter_binding_file(mhc1_files, "filtered_mhc1")
filter_binding_file(mhc2_files, "filtered_mhc2")
