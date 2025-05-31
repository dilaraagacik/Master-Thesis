args <- commandArgs(trailingOnly = TRUE)
input_csv <- args[1]
output_dir <- args[2]

MHC_gt <- read.csv(input_csv, stringsAsFactors = FALSE)

# Format MHC-I
for(a in c("A","B","C")){
  for(i in 1:2) MHC_gt[[paste0(a,".", i)]] <- paste0("HLA-", a, MHC_gt[[paste0(a,".", i)]])
}

# Remove colons from MHC-II
for(i in grep("D", colnames(MHC_gt))) {
  MHC_gt[[i]] <- gsub(":", "", MHC_gt[[i]])
}

# Collect all MHC-I and MHC-II alleles
mhc1_all <- unique(unlist(MHC_gt[, c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2")]))
mhc2_all <- unique(unlist(lapply(1:nrow(MHC_gt), function(i) c(
  paste0("HLA-DPA1", MHC_gt[i, "DPA1.1"], "-DPB1", MHC_gt[i, "DPB1.1"]),
  paste0("HLA-DPA1", MHC_gt[i, "DPA1.1"], "-DPB1", MHC_gt[i, "DPB1.2"]),
  paste0("HLA-DPA1", MHC_gt[i, "DPA1.2"], "-DPB1", MHC_gt[i, "DPB1.1"]),
  paste0("HLA-DPA1", MHC_gt[i, "DPA1.2"], "-DPB1", MHC_gt[i, "DPB1.2"]),
  paste0("HLA-DQA1", MHC_gt[i, "DQA1.1"], "-DQB1", MHC_gt[i, "DQB1.1"]),
  paste0("HLA-DQA1", MHC_gt[i, "DQA1.1"], "-DQB1", MHC_gt[i, "DQB1.2"]),
  paste0("HLA-DQA1", MHC_gt[i, "DQA1.2"], "-DQB1", MHC_gt[i, "DQB1.1"]),
  paste0("HLA-DQA1", MHC_gt[i, "DQA1.2"], "-DQB1", MHC_gt[i, "DQB1.2"]),
  paste0("DRB1_", MHC_gt[i, "DRB1.1"]),
  paste0("DRB1_", MHC_gt[i, "DRB1.2"])
))))

# Write unique alleles to shared files
writeLines(mhc1_all, file.path(output_dir, "mhc1_alleles.txt"))
writeLines(mhc2_all, file.path(output_dir, "mhc2_alleles.txt"))
