library(readr)
library(dplyr)
library(stringr)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
reference_file <- args[2]
output_file <- args[3]

cat("📥 Input file:", input_file, "\n")
cat("📘 Reference file:", reference_file, "\n")

# Load input CSV
hla_df <- read_csv(input_file, col_types = cols(.default = "c"))

# Load reference
ref_df <- read_csv(reference_file, col_types = cols(
  Gene = col_character(),
  Allele = col_character(),
  P_Value = col_character()
)) %>%
  mutate(
    Gene = str_trim(Gene),
    Allele = str_trim(Allele)
  )

# Function to map allele to P_Value
map_allele <- function(column_name, allele_value) {
  if (is.na(allele_value) || allele_value == "") return(allele_value)

  gene <- str_extract(column_name, "^[A-Z]+")
  allele <- allele_value

  # Try to find match in reference file
  matched <- ref_df %>%
    filter(Gene == gene & str_starts(Allele, allele))

  if (nrow(matched) > 0 && !is.na(matched$P_Value[1]) && matched$P_Value[1] != "") {
    return(str_remove(matched$P_Value[1], "P"))
  } else {
    return(allele_value)  # Keep original if P_Value is missing or empty
  }
}

# Apply mapping across all allele columns (excluding sample_id)
mapped_df <- hla_df %>%
  mutate(across(
    -sample_id,
    ~mapply(map_allele, cur_column(), .),
    .names = "{.col}"
  ))

# Write output
write_csv(mapped_df, output_file)
cat("✅ Output written to:", output_file, "\n")
