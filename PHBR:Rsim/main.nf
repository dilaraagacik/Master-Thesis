#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Ensure output dirs exist
new File(params.output_dir).mkdirs()
new File(params.annotation_dir).mkdirs()



// === PROCESSES ===

process AnnotateVCF {
  tag { sample_id }

  input:
  tuple val(sample_id), path(vcf)

  output:
  tuple val(sample_id), path("*_anno.hg38_multianno.txt"), emit: annotated_file

  script:
  """
  filename=\$(basename ${vcf})

  perl ${params.annovar_path}/table_annovar.pl ${vcf} \\
    ${params.annovar_db} \\
    -buildver hg38 \\
    -out ./\${filename}_anno \\
    -protocol refGene \\
    -operation g \\
    -remove -polish -vcfinput -nastring .
  """
}


process GeneratePeptides {
  tag { sample_id }

  input:
    tuple val(sample_id), path(multianno_txt), path(mapping_file), path(mane_file), path(fasta_file)

  output:
    tuple val(sample_id), path("peptides_${sample_id}.txt")


  script:
  """
  Rscript ${params.process_mutations} ${multianno_txt} peptides_${sample_id}.txt ${sample_id} ${mapping_file} ${mane_file} ${fasta_file}
  """
}

process GeneratePeptidesFromMAF {
  tag { sample_id }

  input:
    tuple val(sample_id), path(maf_file), path(mapping_file), path(mane_file), path(fasta_file)

  output:
    tuple val(sample_id), path("peptides_${sample_id}.txt")

  script:
  """
  Rscript ${params.process_mutations_from_maf} ${maf_file} peptides_${sample_id}.txt ${sample_id} ${mapping_file} ${mane_file} ${fasta_file}
  """
}




process PMapping {
  tag "PMapping - Mapping alleles"
  

  input:
  path(genotypes) // Path to the genotype file (containing multiple samples and their alleles)
  path(reference)  // Path to the reference file with allele-P-value mappings

  output:
  path("filtered_hla.csv"), emit: filtered_hla

  script:
  """
  # Run the R script for allele mapping on all peptides (no sample-specific processing)
  Rscript ${params.p_mapping} "${genotypes}" "${reference}" filtered_hla.csv

  echo "SUCCESS: Completed P-Mapping for all alleles."
  """
}
process PreparePeptideFASTAs {

  input:
  tuple val(tcga_id), path(peptide_table)

  output:
  tuple val(tcga_id), path("mhc1_${tcga_id}.fasta"), path("mhc2_${tcga_id}.fasta")

  script:
  """
  # MHC-I FASTA
  awk -F '\\t' '
    NR==1 {
      for(i=1;i<=NF;i++) {
        if(\$i=="Simple_Mutation") mut_col=i
        if(\$i=="Hugo_Symbol") hugo_col=i
        if(\$i=="Peptide_MHC1") mhc1_col=i
        if(\$i=="All_Mutant9mers") nine_col=i
      }
      next
    }
    {
      if (mhc1_col && \$(mhc1_col) != "" && \$(mhc1_col) != "NA") {
        # Prefer Peptide_MHC1 if present
        printf(">%s_%s\\n%s\\n", \$(mut_col), (hugo_col ? \$(hugo_col) : ""), \$(mhc1_col))
      } else if (nine_col && \$(nine_col) != "" && \$(nine_col) != "NA") {
        n=split(\$(nine_col), arr, ";")
        for(j=1;j<=n;j++) {
          if (arr[j] != "" && arr[j] != "NA")
            printf(">%s_%s_%d\\n%s\\n", \$(mut_col), (hugo_col ? \$(hugo_col) : ""), j, arr[j])
        }
      }
    }
  ' ${peptide_table} > mhc1_${tcga_id}.fasta

  # MHC-II FASTA
  awk -F '\\t' '
    NR==1 {
      for(i=1;i<=NF;i++) {
        if(\$i=="Simple_Mutation") mut_col=i
        if(\$i=="Hugo_Symbol") hugo_col=i
        if(\$i=="Peptide_MHC2") mhc2_col=i
        if(\$i=="All_Mutant15mers") fifteen_col=i
      }
      next
    }
    {
      if (mhc2_col && \$(mhc2_col) != "" && \$(mhc2_col) != "NA") {
        printf(">%s_%s\\n%s\\n", \$(mut_col), (hugo_col ? \$(hugo_col) : ""), \$(mhc2_col))
      } else if (fifteen_col && \$(fifteen_col) != "" && \$(fifteen_col) != "NA") {
        n=split(\$(fifteen_col), arr, ";")
        for(j=1;j<=n;j++) {
          if (arr[j] != "" && arr[j] != "NA")
            printf(">%s_%s_%d\\n%s\\n", \$(mut_col), (hugo_col ? \$(hugo_col) : ""), j, arr[j])
        }
      }
    }
  ' ${peptide_table} > mhc2_${tcga_id}.fasta
  """
}






process RunNetMHCpan {
  tag { tcga_id }
  errorStrategy 'ignore'

  input:
    tuple val(tcga_id), path(mhc1_alleles), path(mhc1_fasta)

  output:
    tuple val(tcga_id), path("mhc1_*.xls"), emit: mhc1_results

  script:
  """
  if [ ! -s "${mhc1_fasta}" ]; then
    echo "⚠️ FASTA file is empty for ${tcga_id}, skipping NetMHCpan."
    touch "mhc1_${tcga_id}_DUMMY.xls"
    exit 0
  fi

  any_valid=false
  while read allele_full; do
    if [[ "\$allele_full" == *NA* ]]; then
      echo "⚠️ Skipping invalid allele: \$allele_full"
      # Still touch file so Nextflow is satisfied
      touch "mhc1_${tcga_id}_\${allele_full}.xls"
    else
      any_valid=true
      clean_allele=\$(echo "\$allele_full" | sed "s/_twice//")
      output="mhc1_${tcga_id}_\${allele_full}.xls"
      echo "🔬 Running NetMHCpan for allele: \$clean_allele"
      tcsh ${params.tool_mhc1} -a "\$clean_allele" -f "${mhc1_fasta}" -inptype 0 -l 8,9,10,11 -BA -xls -xlsfile "\$output"
    fi
  done < ${mhc1_alleles}

  # Safety net: if no valid alleles, emit at least one dummy file
  if [ "\$any_valid" = false ]; then
    echo "⚠️ No valid MHC-I alleles found for ${tcga_id}, generating dummy output."
    touch "mhc1_${tcga_id}_DUMMY.xls"
  fi
  """
}


process RunNetMHCIIpan {
  tag { tcga_id }
  errorStrategy 'ignore'

  input:
    tuple val(tcga_id), path(mhc2_alleles), path(mhc2_fasta)

  output:
    tuple val(tcga_id), path("mhc2_*.xls"), emit: mhc2_results

  script:
  """
  if [ ! -s "${mhc2_fasta}" ]; then
    echo "⚠️ FASTA file is empty for ${tcga_id}, skipping NetMHCIIpan."
    touch "mhc2_${tcga_id}_DUMMY.xls"
    exit 0
  fi

  any_valid=false
  while read allele_full; do
    if [[ "\$allele_full" == *NA* ]]; then
      echo "⚠️ Skipping invalid allele: \$allele_full"
      touch "mhc2_${tcga_id}_\${allele_full}.xls"
    else
      any_valid=true
      clean_allele=\$(echo "\$allele_full" | sed "s/_twice//")
      output="mhc2_${tcga_id}_\${allele_full}.xls"
      echo "🔬 Running NetMHCIIpan for allele: \$clean_allele"
      tcsh ${params.tool_mhc2} -a "\$clean_allele" -f "${mhc2_fasta}" -inptype 0 -length 15 -xls -xlsfile "\$output"
    fi
  done < ${mhc2_alleles}

  if [ "\$any_valid" = false ]; then
    echo "⚠️ No valid MHC-II alleles found for ${tcga_id}, generating dummy output."
    touch "mhc2_${tcga_id}_DUMMY.xls"
  fi
  """
}


process CalculatePHBR_MHC1 {
  errorStrategy 'ignore'
  tag { tcga_id }
    publishDir "/home/dilara/projects/perimeta/results/nextflow_pipeline/phbr_nextflow/results/phbr", mode: 'copy'

  input:
    tuple val(tcga_id), path(mhc1_files)

  output:
    path("phbr_${tcga_id}_mhc1.csv")

  script:
  """
  mkdir -p mhc1_files
  valid_count=0

  for f in ${mhc1_files}; do
    if [[ "\$f" != *DUMMY.xls && "\$f" != *NA*.xls ]]; then
      mv "\$f" mhc1_files/
      valid_count=\$((valid_count + 1))
    else
      echo "⚠️ Skipping dummy/NA file: \$f"
    fi
  done

  output_file="phbr_${tcga_id}_mhc1.csv"

  if [ \$valid_count -eq 0 ]; then
    echo "⚠️ No valid MHC-I files found for ${tcga_id}, writing fallback PHBR result."
    echo "sample_id,PHBR_MHC1" > "\$output_file"
    echo "${tcga_id},NA" >> "\$output_file"
  else
    python3 ${params.calculate_phbr} mhc1_files dummy_dir --sample_id ${tcga_id}

    # 🔐 Safety net: ensure output file exists
    if [ ! -f "\$output_file" ]; then
      echo "⚠️ PHBR file not created by script. Writing fallback."
      echo "sample_id,PHBR_MHC1" > "\$output_file"
      echo "${tcga_id},NA" >> "\$output_file"
    fi
  fi
  """

}
process CalculatePHBR_MHC2 {
  tag { tcga_id }
  errorStrategy 'ignore'
  publishDir "/home/dilara/projects/perimeta/results/nextflow_pipeline/phbr_nextflow/results/phbr", mode: 'copy'

  input:
    tuple val(tcga_id), path(mhc2_files)

  output:
    path("phbr_${tcga_id}_mhc2.csv")

  script:
  """
  mkdir -p mhc2_files
  valid_count=0

  for f in ${mhc2_files}; do
    if [[ "\$f" != *DUMMY.xls && "\$f" != *NA*.xls ]]; then
      mv "\$f" mhc2_files/
      valid_count=\$((valid_count + 1))
    else
      echo "⚠️ Skipping dummy/NA file: \$f"
    fi
  done

  output_file="phbr_${tcga_id}_mhc2.csv"

  if [ \$valid_count -eq 0 ]; then
    echo "⚠️ No valid MHC-II files found for ${tcga_id}, writing fallback PHBR result."
    echo "sample_id,PHBR_MHC2" > "\$output_file"
    echo "${tcga_id},NA" >> "\$output_file"
  else
    python3 ${params.calculate_phbr} dummy_dir mhc2_files --sample_id ${tcga_id}

    # Ensure output file is written
    if [ ! -f "\$output_file" ]; then
      echo "⚠️ PHBR script completed but did not produce output — writing fallback file."
      echo "sample_id,PHBR_MHC2" > "\$output_file"
      echo "${tcga_id},NA" >> "\$output_file"
    fi
  fi
  """
}





// ✅ **Define Input Samples (FASTQ Files)**
def read_pairs = nextflow.Channel
    .fromPath("${params.fastq_path}/*/", type: 'dir')
    .filter { it.exists() && it.isDirectory() }
    .map { sampleDir ->
        def sample = sampleDir.getName()

        // Try common FASTQ naming patterns
        def patterns = [
            "${sample}_1.fastq.gz"      : "${sample}_2.fastq.gz",
            "${sample}_R1.fastq.gz"     : "${sample}_R2.fastq.gz",
            "${sample}_read1.fastq.gz"  : "${sample}_read2.fastq.gz",
            "${sample}_R1_001.fastq.gz" : "${sample}_R2_001.fastq.gz"
        ]

        def pair = patterns.find { new File(sampleDir, it.key).exists() && new File(sampleDir, it.value).exists() }

        if (!pair) {
            log.warn "⚠️ No matching FASTQ pair found for sample: ${sample} in ${sampleDir}"
            return null
        }

        def read1 = file("${sampleDir}/${pair.key}")
        def read2 = file("${sampleDir}/${pair.value}")
        return tuple(sample, read1, read2)
    }
    .filter { it != null }


process BamConversion {
    tag "BAM ${sample_id}"

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz"), emit: paired_fastq

    script:
    """
    ${params.bam_to_fastq_script} "${sample_id}" "${bam_file}" ${params.threads}
    """
}
process HLAHD {
    tag "HLA-HD - ${sample_id}"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.HLAHD.results.txt"), emit: hla_genotype

    script:
    """
    echo "🔬 Running HLA-HD for sample: ${sample_id}"

    source ${params.conda_profile}
    conda activate "${params.hlahd_env}"

    toolfolder="${params.hlahd_toolfolder}"
    toolbin="\$toolfolder/bin"
    export PATH="\$toolbin:\$PATH"

    hlahd.sh -t ${params.threads} -m 30 -c 0.95 -f "\$toolfolder/freq_data/" \
      "${read1}" "${read2}" "\$toolfolder/HLA_gene.split.txt" "\$toolfolder/dictionary/" \
      "${sample_id}" "."

    result_file=\$(find . -type f -path "*/result/*.txt" | head -n 1)

    if [[ -s "\$result_file" ]]; then
        cp "\$result_file" "${sample_id}.HLAHD.results.txt"
    else
        echo "❌ No result file found for sample: ${sample_id}" >&2
        exit 1
    fi
    """
}

process ParseHLAHDToCSV {
    tag "Parse HLA-HD Genotypes"

    input:
    path hla_txts

    output:
    path("hla_genotype_table.csv"), emit: parsed_genotypes

    script:
    """
    bash ${params.parse_hlahd_to_csv_script}
    """
}


process ArcasHLA_FASTQ {
    tag "ArcasHLA_FASTQ - ${sample_id}"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_1.genotype.json"), emit: hla_genotype_fastq

    script:
    """
    echo "🚀 Running ArcasHLA (FASTQ) on ${sample_id}"
    "${params.toolscript_arcasHLA}" genotype "${read1}" "${read2}" \
      -g A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1 -o "." -t ${params.threads} -v
    """
}

process ArcasHLA_BAM {
  tag "ArcasHLA_BAM - ${sample_id}"

  input:
    tuple val(sample_id), path(bam_file)

  output:
    tuple val(sample_id), path("${sample_id}_1.genotype.json"), emit: hla_genotype_bam

  script:
  """
  echo "🔄 Sorting BAM by name for ${sample_id}"
  samtools sort -n -@ ${params.threads} -o ${sample_id}_sorted.bam ${bam_file}

  echo "🔄 Converting BAM to FASTQ for ${sample_id}"
  samtools fastq -@ ${params.threads} \\
    -1 ${sample_id}_1.fastq.gz \\
    -2 ${sample_id}_2.fastq.gz \\
    -0 /dev/null -s /dev/null -n ${sample_id}_sorted.bam

  echo "🚀 Running ArcasHLA (BAM) on ${sample_id}"
  "${params.toolscript_arcasHLA}" genotype "${sample_id}_1.fastq.gz" "${sample_id}_2.fastq.gz" \\
    -g A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1 -o "." -t ${params.threads} -v
  """
}



process PMappingAndFilter {
    tag "PMappingAndFilter - ${sample_id}"

    input:
    tuple val(sample_id), path(hla_genotype)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_genotype.json"), emit: filtered_hla_genotype

    script:
    """
    # Load Conda environment with R

    # Run the R script
    Rscript ${params.tool_pmapping} "${hla_genotype}" "${params.p_mapping_dataset}" "${sample_id}_filtered_genotype.json"

    echo "✅ SUCCESS: Completed P-Mapping for ${sample_id}"
    """
}

process ParseGenotypesToCSV {
  tag "Parse Genotypes"

  input:
  tuple val(sample_id), path(json_file)


  output:
  path("filtered_hla.csv"), emit: parsed_genotypes

  script:
  """
  echo "sample_id,A.1,A.2,B.1,B.2,C.1,C.2,DPA1.1,DPA1.2,DPB1.1,DPB1.2,DQA1.1,DQA1.2,DQB1.1,DQB1.2,DRB1.1,DRB1.2" > filtered_hla.csv

  for file in *.json; do
    sample_id=\$(basename \$file | sed 's/_filtered_genotype.json//')

    alleles=\$(jq -r '
      [
        .A[0], .A[1],
        .B[0], .B[1],
        .C[0], .C[1],
        .DPA1[0], .DPA1[1],
        .DPB1[0], .DPB1[1],
        .DQA1[0], .DQA1[1],
        .DQB1[0], .DQB1[1],
        .DRB1[0], .DRB1[1]
      ]
      | map(sub(".*\\\\*"; "") | split(":")[0:2] | join(":"))
      | join(",")
    ' \$file)

    echo "\$sample_id,\$alleles" >> filtered_hla.csv
  done
  """
}



process SplitHLAByPatient {
  tag "Split HLA by patient"
  input:
    path "filtered_hla.csv"

  output:
    path("mhc1_*.txt"), emit: mhc1
    path("mhc2_*.txt"), emit: mhc2

  script:
  """
  tail -n +2 filtered_hla.csv | cut -d',' -f1 | cut -d'-' -f1-3 | while read short_id; do
    echo "🔵 Processing short ID: \$short_id"
    clean_id=\$(echo "\$short_id" | tr -d '"')
    awk -F',' -v id="\$short_id" 'NR==1 || \$1 == id' filtered_hla.csv > \${short_id}_filtered.csv

    if [ -s \${short_id}_filtered.csv ]; then
      echo "✅ Found data for \$short_id"
      Rscript ${params.convert_hla} \${short_id}_filtered.csv .

      if [[ -s mhc1_alleles.txt ]]; then
        awk '\$1 != "NA" { count[\$1]++ } 
             END {
               for (a in count) {
                 if (count[a]>1)
                   print a"_twice"
                 else
                   print a
               }
             }' mhc1_alleles.txt > mhc1_\${clean_id}.txt
      fi

      if [[ -s mhc2_alleles.txt ]]; then
        awk '\$1 != "NA" { count[\$1]++ } 
             END {
               for (a in count) {
                 if (count[a]>1)
                   print a"_twice"
                 else
                   print a
               }
             }' mhc2_alleles.txt > mhc2_\${clean_id}.txt
      fi

      rm \${short_id}_filtered.csv
      rm mhc2_alleles.txt
      rm mhc1_alleles.txt
    else
      echo "⚠️ No matching rows found for \$short_id"
    fi
  done
  """
}



process detect_bam_type {
  tag { real_sample_id }

  input:
    tuple val(sample_id), val(real_sample_id), path(bam_file)

  output:
    tuple val(sample_id), val(real_sample_id), path(bam_file), val(is_paired)

  script:
  """
  if samtools view -f 1 ${bam_file} | head -n 1 | grep -q .; then
    echo "paired"
  else
    echo "single"
  fi
  """
}


process featureCounts_quantification {
  tag { real_sample_id }

  input:
    tuple val(sample_id), val(real_sample_id), path(bam_file)

  output:
    path("${real_sample_id}_counts.txt")

  script:
  """
  echo "🔍 Inspecting BAM for paired-end layout: ${real_sample_id}"
  if samtools view -f 1 ${bam_file} | head -n 1 | grep -q .; then
    echo "📦 Running featureCounts in PAIRED mode for ${real_sample_id}"
    featureCounts \\
      -T ${params.threads} \\
      -p \\
      -M \\
      --primary \\
      -t exon \\
      -g gene_id \\
      -a ${params.gtf} \\
      -o ${real_sample_id}_counts.txt \\
      ${bam_file}
  else
    echo "📦 Running featureCounts in SINGLE mode for ${real_sample_id}"
    featureCounts \\
      -T ${params.threads} \\
      -M \\
      --primary \\
      -t exon \\
      -g gene_id \\
      -a ${params.gtf} \\
      -o ${real_sample_id}_counts.txt \\
      ${bam_file}
  fi

  echo "✅ Done: ${real_sample_id}"
  cat ${real_sample_id}_counts.txt.summary || echo "⚠️ No summary file found."
  """
}



process merge_and_normalize_counts {
  publishDir "results/featureCounts", mode: 'copy'
    tag "Merging and normalizing featureCounts outputs"

    input:
    path(counts_files)

    output:
    path("merged_counts.tsv"), emit: merged_counts_file


    script:
    """
    counts_dir=\$(dirname \$(ls ${counts_files} | head -n1))
    echo "Using counts directory: \$counts_dir"

    Rscript ${params.rscript_merge_norm} \\
      \$counts_dir merged_counts.tsv
    """
}


process detect_zero_expression_genes {

    input:
    path merged_counts_file

    output:
    path "zero_genes_per_patient/*.txt"

    script:
    """
    mkdir zero_genes_per_patient
    Rscript ${params.rscript_zero_genes} ${merged_counts_file} zero_genes_per_patient ${params.biomart_csv}
    """
}

process filter_binding_by_expression {
    tag { sample_id }

    input:
    tuple val(sample_id), path(mhc1_files), path(mhc2_files), path(zero_genes_file)

    output:
    tuple val(sample_id), path("filtered_mhc1_*.xls"), path("filtered_mhc2_*.xls")

    script:
    """
    Rscript ${params.rscript_filter_binding} \\
      --mhc1_files ${mhc1_files.join(',')} \\
      --mhc2_files ${mhc2_files.join(',')} \\
      --zero_genes ${zero_genes_file} \\
      --driver_genes ${params.driver_genes} \\
      --out_prefix ${sample_id}
    """
}

process CalculateRexpSample {
  tag { sample_id }
  errorStrategy 'ignore'

  input:
    tuple val(sample_id), path(mhc1_file), path(mhc2_file)

  output:
    path "rexp_${sample_id}.csv"

  script:
script:
"""
Rscript ${params.r_script_rexp_sample} \\
  --sample_id ${sample_id} \\
  --mhc1_file ${mhc1_file} \\
  --mhc2_file ${mhc2_file} \\
  --signature_file ${params.signature_file} \\
  --core_rdata ${params.rand_core_rdata} \\
  --metadata_file ${params.metadata} \\
  --output rexp_${sample_id}.csv
"""

}

process MergeRexpResults {
  tag "Merging Rexp Results"

  input:
    path "rexp_*.csv"

  output:
    path "rexp_output.csv", emit: rexp_out

  script:
  """
  # 🛠️ Merge all rexp CSVs
  first=\$(ls rexp_*.csv | head -n 1)
  head -n 1 \$first > rexp_output_temp.csv
  for f in rexp_*.csv; do tail -n +2 \$f >> rexp_output_temp.csv; done

  # 🧹 Clean duplicated rows and remove unwanted 'alleles' rows
  awk -F',' 'NR==1 || (!seen[\$1]++ && \$1 != "alleles")' rexp_output_temp.csv > rexp_output.csv

  echo "✅ Cleaned and merged rexp_output.csv generated."
  """
}


process CalculateRobsSample {
  cache = false
  tag { sample_id }


  input:
    tuple val(sample_id), val(mhc1_files), val(mhc2_files)

  output:
    path("robs_${sample_id}.csv")

  script:
    def mhc1_str = mhc1_files.join(',')
    def mhc2_str = mhc2_files.join(',')

    """
    echo "📥 Sample: ${sample_id}"
    echo "🔗 MHC1 files: ${mhc1_str}"
    echo "🔗 MHC2 files: ${mhc2_str}"

    Rscript ${params.r_script_robs} \\
      --mhc1_files "${mhc1_str}" \\
      --mhc2_files "${mhc2_str}" \\
      --output robs_${sample_id}.csv
    """
}

process MergeRobsResults {
  input:
    path "robs_*.csv"

  output:
    path "binding_summary.csv", emit: robs_out

  script:
  """
  first=\$(ls robs_*.csv | head -n 1)
  head -n 1 "\$first" > binding_summary.csv
  for f in robs_*.csv; do tail -n +2 "\$f" >> binding_summary.csv; done
  """
}

process JoinRobsRexp {
  publishDir "/home/dilara/projects/perimeta/results/nextflow_pipeline/phbr_nextflow/results/Rsim", mode: 'copy'

  input:
    path robs_csv
    path rexp_csv

  output:
    path "Robs_Rexp_Rratio_summary.csv", emit: rsim_summary

  script:
  """
  Rscript -e '
    library(readr)
    library(dplyr)

    rob <- read_csv("binding_summary.csv") %>%
      distinct(sample_id, .keep_all = TRUE)

    rexp <- read_csv("rexp_output.csv") %>%
      distinct(sample_id, .keep_all = TRUE)

    final <- left_join(rob, rexp, by = "sample_id") %>%
      mutate(
        Rratio_MHC1 = Robs_MHC1 / Rexp_MHC1,
        Rratio_MHC2 = Robs_MHC2 / Rexp_MHC2
      )

    write_csv(final, "Robs_Rexp_Rratio_summary.csv")
  '
  """
}



workflow {

    // === Define VCF input files if input_format is "vcf"
  def vcf_files
  if (params.input_format == 'vcf') {
    vcf_files = Channel
      .fromPath("${params.vcf_dir}/**/*.vcf.gz", type: 'file')
      .map { file ->
        def sample_id = file.getBaseName().tokenize('.')[0]
        tuple(sample_id, file)
      }
  }


    // === Load file-to-barcode mapping from CSV ===
  def maf_files

  if (params.barcode_map_csv) {
      log.info "📋 Reading MAF file paths from barcode map CSV: ${params.barcode_map_csv}"

      def barcode_map_list = []

      new File(params.barcode_map_csv).splitEachLine(",") { line ->
          if (line.size() >= 3 && !line[0].toLowerCase().contains("file")) {
              try {
                  def file_path = line[0].replaceAll('"', '')
                  def barcode = line[2].replaceAll('"', '').tokenize('-')[0..2].join('-')
                  barcode_map_list << tuple(barcode, file(file_path))
              } catch (Exception e) {
                  log.warn "⚠️ Skipping malformed line: ${line} -- ${e.message}"
              }
          }
      }

      maf_files = nextflow.Channel.fromList(barcode_map_list)




        
  } else {
      log.info "📁 No barcode_map_csv provided — using files from maf_dir: ${params.maf_dir}"

      maf_files = nextflow.Channel
          .fromPath("${params.maf_dir}/*/*.maf.gz", type: 'file')
          .map { file ->
              def sample_id = file.getBaseName().tokenize('.')[0]
              tuple(sample_id, file)
          }
      
  }

  if (!["phbr", "rsim"].contains(params.mode)) {
    error "Invalid mode: '${params.mode}'. Please set --mode phbr or --mode rsim."
  }

  if (!["vcf", "maf"].contains(params.input_format)) {
    error "Invalid input format: '${params.input_format}'. Please set --input_format vcf or maf."
  }

  // === Step 1: Annotate or Directly Use MAF and Generate Peptides

  def peptides
  if (params.input_format == "vcf") {
    def annotated = AnnotateVCF(vcf_files)
    peptides = annotated.map { sample_id, anno_file ->
      def id = (anno_file.getBaseName() =~ /([a-f0-9\-]{36})/)[0][1]
      tuple(id, anno_file, file(params.mapping_file), file(params.mane_file), file(params.fasta_file))
    } | GeneratePeptides

  } else {
    peptides = maf_files.map { sample_id, maf ->
      tuple(sample_id, maf, file(params.mapping_file), file(params.mane_file), file(params.fasta_file))
    } | GeneratePeptidesFromMAF
  }


  def fastas = PreparePeptideFASTAs(peptides)
  def mhc1_fasta_channel = fastas.map { id, fasta1, fasta2 -> tuple(id, fasta1) }
  def mhc2_fasta_channel = fastas.map { id, fasta1, fasta2 -> tuple(id, fasta2) }

  def hla_split_results
  def filtered_hla_csv

  if (params.genotypes) {
    log.info "📥 Using provided filtered HLA genotype CSV"
    filtered_hla_csv = file(params.genotypes)

  } else if (params.data_path) {
    log.info "🧬 No genotype CSV provided — running HLA caller"

    if (params.hla_caller == 'arcasHLA') {
      def arcas_results

      if (params.hla_input_type == 'fastq') {
        if (!params.fastq_path) error "❌ Missing --fastq_path for FASTQ input"
        arcas_results = read_pairs | ArcasHLA_FASTQ

      } else if (params.hla_input_type == 'bam') {
        def bam_files = nextflow.Channel
          .fromPath("${params.data_path}/**/*.bam", type: 'file')
          .map { file -> tuple(file.getBaseName().tokenize('.')[0], file) }

        arcas_results = bam_files | ArcasHLA_BAM

      } else {
        error "❌ Invalid --hla_input_type '${params.hla_input_type}'. Use 'fastq' or 'bam'."
      }

      def pmapped = arcas_results | PMappingAndFilter
      def parsed = pmapped.map { sid, json -> tuple(sid, json) } | ParseGenotypesToCSV
      filtered_hla_csv = parsed.parsed_genotypes
    }

    else if (params.hla_caller == 'hlahd') {
      def hlahd_inputs

      if (params.hla_input_type == 'fastq') {
        if (!params.fastq_path) error "❌ Missing --fastq_path for FASTQ input"
        hlahd_inputs = read_pairs

      } else if (params.hla_input_type == 'bam') {
        def bam_files = nextflow.Channel
          .fromPath("${params.data_path}/**/*.bam", type: 'file')
          .map { file -> tuple(file.getBaseName().tokenize('.')[0], file) }

        hlahd_inputs = bam_files | BamConversion
      } else {
        error "❌ Invalid --hla_input_type '${params.hla_input_type}' for HLA-HD. Use 'fastq' or 'bam'."
      }

      def hlahd_results = hlahd_inputs | HLAHD
      filtered_hla_csv = hlahd_results.map { it -> it[1] } | ParseHLAHDToCSV
    }

    else {
      error "❌ Unknown --hla_caller '${params.hla_caller}'. Use 'arcasHLA' or 'hlahd'."
    }

  } else {
    error "❌ No valid HLA input provided. Use --genotypes or --data_path."
  }




  hla_split_results = SplitHLAByPatient(filtered_hla_csv)

  def hla_alleles_mhc1 = hla_split_results.mhc1.flatten().map {
    file -> tuple(file.getBaseName().replace("mhc1_", "").replace(".txt", ""), file)
  }
  def hla_alleles_mhc2 = hla_split_results.mhc2.flatten().map {
    file -> tuple(file.getBaseName().replace("mhc2_", "").replace(".txt", ""), file)
  }

  def netmhc1 = hla_alleles_mhc1.join(mhc1_fasta_channel)
    .map { id, alleles, fasta -> tuple(id, alleles, fasta) } | RunNetMHCpan

  def netmhc2 = hla_alleles_mhc2.join(mhc2_fasta_channel)
    .map { id, alleles, fasta -> tuple(id, alleles, fasta) } | RunNetMHCIIpan

  // === PHBR Path ===
  if (params.mode == "phbr") {
    netmhc1.groupTuple().map { id, files -> tuple(id, files.flatten()) } | CalculatePHBR_MHC1
    netmhc2.groupTuple().map { id, files -> tuple(id, files.flatten()) } | CalculatePHBR_MHC2
  }

  // === RSIM Path ===
  if (params.mode == "rsim") {

    def merged_counts

    if (params.merged_counts) {
      log.info "📊 Using provided merged counts: ${params.merged_counts}"
      merged_counts = nextflow.Channel.fromPath(params.merged_counts)
    } else {
      def bam_files_rsim = nextflow.Channel
        .fromPath("${params.data_path}/**/*.bam")
        .map { bam -> tuple(bam.getBaseName(), bam.getParent().getName(), bam) }

      def counts = bam_files_rsim | featureCounts_quantification
      merged_counts = counts.collect() | merge_and_normalize_counts
    }

    def zero_genes_lists = merged_counts | detect_zero_expression_genes
    def zero_genes_mapped = zero_genes_lists.flatten().map { file ->
      def sample_id = file.getBaseName().replace("_zero_genes", "").replace(".txt", "")
      tuple(sample_id, file)
    }

    def mhc1_grouped = netmhc1.groupTuple()
    def mhc2_grouped = netmhc2.groupTuple()

    def filtered_mhc = mhc1_grouped
      .join(mhc2_grouped, by: 0)
      .join(zero_genes_mapped, by: 0)
      .map { id, mhc1_files, mhc2_files, zero_genes_file ->
        tuple(id, mhc1_files.flatten(), mhc2_files.flatten(), zero_genes_file)
      } | filter_binding_by_expression

    def filtered_robs_inputs = filtered_mhc.map { id, mhc1, mhc2 -> tuple(id, [mhc1], [mhc2]) }

    def allele_files = hla_split_results.mhc1.flatten()
      .map { f -> tuple(f.getBaseName().replace("mhc1_", "").replace(".txt", ""), f) }
      .join(hla_split_results.mhc2.flatten().map { f -> tuple(f.getBaseName().replace("mhc2_", "").replace(".txt", ""), f) }, by: 0)

    def rexp_inputs = allele_files.map { sample_id, mhc1_file, mhc2_file ->
      tuple(sample_id, mhc1_file, mhc2_file)
    }

    def rexp_results = rexp_inputs | CalculateRexpSample
    def rexp_merged = rexp_results.collect() | MergeRexpResults

    def robs_results = filtered_robs_inputs | CalculateRobsSample
    def robs_merged = robs_results.collect() | MergeRobsResults

    JoinRobsRexp(robs_merged.robs_out, rexp_merged.rexp_out)
  }
}
