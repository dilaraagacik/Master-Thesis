# Load libraries
library(jsonlite)
library(dplyr)
library(readr)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
hla_genotype_file <- args[1]
reference_file <- args[2]
output_file <- args[3]

cat("📂 Loading HLA genotype JSON from:", hla_genotype_file, "\n")
cat("📂 Loading reference CSV from:", reference_file, "\n")

# Load JSON HLA data
hla_data <- tryCatch({
  fromJSON(hla_genotype_file)
}, error = function(e) {
  stop("❌ ERROR: Failed to read JSON file - ", e$message)
})

cat("✅ Loaded JSON file successfully. Number of genes:", length(hla_data), "\n")

# Load the reference dataset
reference_expanded <- tryCatch({
  read_csv(reference_file, col_types = cols(
    Gene = col_character(),
    Allele = col_character(),
    P_Value = col_character()
  )) %>%
    mutate(Gene = trimws(Gene), Allele = trimws(Allele))
}, error = function(e) {
  stop("❌ ERROR: Failed to read reference CSV - ", e$message)
})

cat("✅ Loaded reference file successfully. Number of records:", nrow(reference_expanded), "\n")

# Function to map allele to P_Value
map_allele_to_pvalue <- function(allele) {
  parts <- strsplit(allele, "\\*")[[1]]
  if (length(parts) < 2) return(allele)  # Return original if malformed

  gene <- parts[1]
  allele_value <- parts[2]

  match <- reference_expanded %>%
    filter(Gene == gene & Allele == allele_value)

  if (nrow(match) > 0) {
    p_value <- gsub("P", "", match$P_Value)
    return(paste(gene, p_value, sep = "*"))
  } else {
    return(allele)
  }
}

# Apply mapping function to all alleles
filtered_hla_data <- lapply(hla_data, function(alleles) {
  sapply(alleles, map_allele_to_pvalue)
})

cat("✅ Completed allele mapping. Writing output JSON...\n")

# Write JSON output
write_json(filtered_hla_data, output_file, pretty = TRUE, auto_unbox = TRUE)

# Confirm file creation
if (file.exists(output_file)) {
  cat("✅ SUCCESS: Filtered HLA genotype saved as", output_file, "\n")
} else {
  stop("❌ ERROR: Output file was not created.")
}
