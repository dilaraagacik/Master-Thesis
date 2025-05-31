#!/usr/bin/env Rscript

# === Libraries ===
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(jsonlite)

# === Input arguments from Nextflow ===
args <- commandArgs(trailingOnly = TRUE)
input_file    <- args[1]
output_file   <- args[2]
sample_id     <- args[3]
mapping_file  <- args[4]
mane_file     <- args[5]
fasta_file    <- args[6]


# === Load RefSeq ↔ Ensembl Mapping ===
mapping <- read.delim(mapping_file, stringsAsFactors = FALSE)
mapping$refseq_mrna <- sub("\\.\\d+$", "", mapping$refseq_mrna)

# === Load MANE Select Canonical Transcript Map ===
header_line <- readLines(mane_file, n = 1)
mane_df <- read.delim(mane_file, header = FALSE, sep = "\t", comment.char = "")
colnames(mane_df) <- strsplit(sub("^#", "", header_line), "\t")[[1]]

mane_df <- subset(mane_df, MANE_status == "MANE Select")
mane_df <- mane_df[, c("symbol", "RefSeq_nuc", "Ensembl_nuc")]
mane_df$RefSeq_nuc <- sub("\\.\\d+$", "", mane_df$RefSeq_nuc)
mane_df$Ensembl_nuc <- sub("\\.\\d+$", "", mane_df$Ensembl_nuc)
colnames(mane_df) <- c("gene", "refseq", "ensembl")

canonical_transcripts <- setNames(mane_df$refseq, mane_df$gene)

# === Load Ensembl protein FASTA ===
ensembl_proteins <- readAAStringSet(fasta_file)
transcript_ids <- sub(".*transcript:([A-Z0-9]+)(\\.\\d+)?\\b.*", "\\1", names(ensembl_proteins))
protein_lookup <- setNames(as.character(ensembl_proteins), transcript_ids)

missing_transcripts <- c()

get_protein_seq <- function(ensembl_transcript_id) {
  ensembl_transcript_id <- sub("\\.\\d+$", "", ensembl_transcript_id)
  if (is.na(ensembl_transcript_id) || ensembl_transcript_id == "") return(NULL)
  if (!(ensembl_transcript_id %in% names(protein_lookup))) {
    missing_transcripts <<- unique(c(missing_transcripts, ensembl_transcript_id))
    return(NULL)
  }
  return(protein_lookup[[ensembl_transcript_id]])
}

apply_aa_change <- function(seq, pos, new_aa) {
  if (pos > nchar(seq)) return(NULL)
  substr(seq, pos, pos) <- new_aa
  return(seq)
}

generate_peptides <- function(mutated_seq, position, to_aa, from_aa, gene, mut_id) {
  n <- nchar(mutated_seq)
  message("🧬 Generating peptides for ", mut_id, " | Mutation at position: ", position, " | Protein length: ", n)

  # MHC-I peptide (21-mer)
  pep_mhc1 <- if (position > 10 && position + 10 <= n) substr(mutated_seq, position - 10, position + 10) else NA
  if (is.na(pep_mhc1)) {
    message("⚠️ Cannot generate MHC-I peptide for ", mut_id)
  } else {
    message("✅ MHC-I peptide: ", pep_mhc1)
  }

  # MHC-II peptide (29-mer)
  pep_mhc2 <- if (position > 14 && position + 14 <= n) substr(mutated_seq, position - 14, position + 14) else NA
  if (!is.na(pep_mhc2)) {
    center <- 15
    substr(pep_mhc2, center, center) <- to_aa  # Always apply the mutation at center
    message("✅ MHC-II peptide: ", pep_mhc2)
  } else {
    message("⚠️ Cannot generate MHC-II peptide for ", mut_id)
  }

  return(list(Peptide_MHC1 = pep_mhc1, Peptide_MHC2 = pep_mhc2))
}

process_mutations <- function(mutation, mapping) {
  aa_changes <- strsplit(mutation, ",")[[1]]
  all_peptides <- list()

  for (aa_change in aa_changes) {
    parts <- strsplit(aa_change, ":")[[1]]
    if (length(parts) >= 5) {
      refseq_id <- sub("\\.\\d+$", "", parts[2])
      gene <- parts[1]

      canonical_id <- canonical_transcripts[gene]

      if (is.na(canonical_id)) {
        message("\u26A0\uFE0F Gene ", gene, " not found in MANE Select — skipping")
        next
      }

      if (canonical_id != refseq_id) {
        message("\u23ED Skipping non-canonical transcript: ", gene, " → ", refseq_id)
        next
      }

      aa_mutation <- parts[5]
      match <- regexec("p\\.([A-Z*]+)(\\d+)([A-Z*]+)", aa_mutation)
      match <- regmatches(aa_mutation, match)

      if (length(match[[1]]) > 0) {
        from_aa <- match[[1]][2]
        position <- as.integer(match[[1]][3])
        to_aa <- match[[1]][4]
        if (from_aa == to_aa) next

        ensembl_id <- mapping$ensembl_transcript_id[mapping$refseq_mrna == refseq_id]
        if (length(ensembl_id) == 0) {
          message("❌ No Ensembl transcript for RefSeq ", refseq_id)
          next
        }

        protein_seq <- get_protein_seq(ensembl_id[1])
        if (is.null(protein_seq) || position > nchar(protein_seq)) next

        message("🔍 Mutation: ", gene, " ", aa_mutation, " at pos ", position, " | Protein length: ", nchar(protein_seq))

        original_aa <- substr(protein_seq, position, position)
        message("🔬 Checking AA: expected ", from_aa, ", found ", original_aa)

        if (is.na(original_aa) || original_aa != from_aa) {
          message("⚠️ AA mismatch, skipping: ", gene, " ", aa_mutation)
          next
        }

        mutated_seq <- apply_aa_change(protein_seq, position, to_aa)
        if (!is.null(mutated_seq)) {
          mut_id <- paste0(gene, "_", sub("p\\.", "", aa_mutation))
          peptides <- generate_peptides(mutated_seq, position, to_aa, from_aa, gene, mut_id)
          all_peptides <- append(all_peptides, list(data.frame(
            Gene = gene,
            RefSeq = refseq_id,
            AAChange = aa_mutation,
            Simple_Mutation = mut_id,
            Peptide_MHC1 = peptides$Peptide_MHC1,
            Peptide_MHC2 = peptides$Peptide_MHC2
          )))
        }
      }
    }
  }

  return(all_peptides)
}

process_anno_file <- function(file_path, mapping) {
  df <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df_missense <- df %>% filter(!is.na(AAChange.refGene), AAChange.refGene != ".")

  unique_mutations <- unique(df_missense$AAChange.refGene)
  all_peptides <- list()
  for (mutation in unique_mutations) {
    mutation_peptides <- process_mutations(mutation, mapping)
    all_peptides <- c(all_peptides, mutation_peptides)
  }

  if (length(all_peptides) > 0) {
    message("✅ Finished processing file: ", file_path, " | Total peptides: ", length(all_peptides))
    return(bind_rows(all_peptides))
  } else {
    message("⚠️ No peptides generated for file: ", file_path)
    return(data.frame())
  }
}

# === Run Main ===
message("🔍 Processing file for sample: ", sample_id)
peptides_df <- process_anno_file(input_file, mapping)

if (nrow(peptides_df) > 0) {
  peptides_df <- peptides_df %>%
    filter(!is.na(Peptide_MHC1) | !is.na(Peptide_MHC2)) %>%
    distinct(Peptide_MHC1, Peptide_MHC2, .keep_all = TRUE)

  peptides_df$sample_id <- sample_id

  base_input <- tools::file_path_sans_ext(basename(input_file))
  output_name <- paste0("peptides_", sample_id, "_", base_input, ".txt")

  write.table(peptides_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("✅ Peptides written to: ", output_name)
} else {
  message("⚠️ No valid peptides found for sample: ", sample_id)
}

if (length(missing_transcripts) > 0) {
  writeLines(unique(missing_transcripts), "missing_transcripts.txt")
}
