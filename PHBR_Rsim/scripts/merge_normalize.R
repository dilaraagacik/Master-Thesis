# merge_and_normalize_safe_final.R

args <- commandArgs(trailingOnly=TRUE)
input_folder <- args[1]
output_counts <- args[2]
output_cpm <- args[3]

library(dplyr)
library(edgeR)

# List all *_counts.txt files
files <- list.files(path = input_folder, pattern = "_counts.txt$", full.names = TRUE)

all_counts <- list()

for (f in files) {
    df <- read.delim(f, comment.char = "#")

    if (!"Geneid" %in% colnames(df)) {
        warning(paste("Skipping file (no Geneid column):", f))
        next
    }

    sample_name <- gsub("_counts.txt$", "", basename(f))

    counts_col <- colnames(df)[ncol(df)]  # Always take LAST column
    df_small <- df[, c("Geneid", counts_col)]
    colnames(df_small)[2] <- sample_name

    all_counts[[sample_name]] <- df_small
}

# Check if any valid files were found
if (length(all_counts) == 0) {
    stop("No valid counts files found. Exiting.")
}

# Merge all
merged_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), all_counts)
merged_counts[is.na(merged_counts)] <- 0

# Save raw counts
write.table(merged_counts, file = output_counts, sep = "\t", quote = FALSE, row.names = FALSE)

