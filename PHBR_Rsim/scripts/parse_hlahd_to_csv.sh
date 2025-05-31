#!/usr/bin/env bash
set -euo pipefail

echo "sample_id,A.1,A.2,B.1,B.2,C.1,C.2,DPA1.1,DPA1.2,DPB1.1,DPB1.2,DQA1.1,DQA1.2,DQB1.1,DQB1.2,DRB1.1,DRB1.2" > hla_genotype_table.csv

for file in *.txt; do
  echo "📂 Processing file: $file"
  sample_id=$(basename "$file" | sed 's/.HLAHD.results.txt//')
  echo "🧬 Sample ID: $sample_id"

  declare -A alleles
  for gene in A B C DPA1 DPB1 DQA1 DQB1 DRB1; do
    alleles[$gene.1]=""
    alleles[$gene.2]=""
  done

  while read -r line; do
    allele1=$(echo "$line" | cut -f1)
    allele2=$(echo "$line" | cut -f2)

    gene=$(echo "$allele1" | cut -d',' -f1 | cut -d'*' -f1 | cut -d'-' -f2)

    echo "  → Gene: $gene | Allele1: $allele1 | Allele2: $allele2"

    if [[ -n "$gene" ]]; then
    alleles[$gene.1]=$(echo "$allele1" | cut -d',' -f1 | sed 's/^.*\*//' | cut -d':' -f1,2)
    alleles[$gene.2]=$(echo "$allele2" | sed 's/^.*\*//' | cut -d':' -f1,2)
    fi

  done < <(grep -Ev '^#' "$file")  # ✅ preserves variable scope

  row="$sample_id"
  for gene in A B C DPA1 DPB1 DQA1 DQB1 DRB1; do
    row+=","${alleles[$gene.1]:-NA}
    row+=","${alleles[$gene.2]:-NA}
  done

  echo "$row" >> hla_genotype_table.csv
done

echo "✅ Finished parsing HLA-HD output to hla_genotype_table.csv"
