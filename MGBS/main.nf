#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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

// ✅ **Step 1: Run ArcasHLA**
process ArcasHLA {
    tag "ArcasHLA - ${sample_id}"
    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_1.genotype.json"), emit: hla_genotype


    script:
    """
    # Log the files being processed
    echo "Processing sample: ${sample_id}, with ${read1} and ${read2}"
    mkdir -p "${params.res_folder}/${sample_id}/"
    # Load Conda environment
    source ${params.conda_profile}
    conda activate ${params.conda_env}

    # Run arcasHLA for genotype
    "${params.toolscript}" genotype "${read1}" "${read2}" -g A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1 -o "." -t ${params.threads} -v
    cp "${sample_id}_1.genotype.json" "${params.res_folder}/${sample_id}/"
    
    # Check the status of arcasHLA
    if [ \$? -eq 0 ]; then
        echo "arcasHLA completed successfully for ${sample_id}"
    else
        echo "❌ Error: arcasHLA failed for ${sample_id}"
    fi
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
    source ${params.conda_profile}
    conda activate ${params.conda_env}

    # Run the R script
    Rscript ${params.tool_pmapping} "${hla_genotype}" "${params.p_mapping_dataset}" "${sample_id}_filtered_genotype.json"

    echo "✅ SUCCESS: Completed P-Mapping for ${sample_id}"
    """
}

process ParseGenotypesToCSV {
    tag "Parse Genotypes"

    input:
    path filtered_jsons

    output:
    path("hla_genotype_table.csv"), emit: parsed_genotypes

    script:
    """
    echo "sample_id,A.1,A.2,B.1,B.2,C.1,C.2,DPA1.1,DPA1.2,DPB1.1,DPB1.2,DQA1.1,DQA1.2,DQB1.1,DQB1.2,DRB1.1,DRB1.2" > hla_genotype_table.csv

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

      echo "\$sample_id,\$alleles" >> hla_genotype_table.csv
    done
    """
}







process ExtractUniqueAlleles {
    tag "Extract Unique MHC Alleles"

    input:
    path genotype_csv
    path r_script

    output:
    path "mhc1_alleles.txt", emit: mhc1_alleles
    path "mhc2_alleles.txt", emit: mhc2_alleles

    script:
    """
    mkdir -p output
    Rscript ${params.convert_hla} ${genotype_csv} output

    mv output/mhc1_alleles.txt ./
    mv output/mhc2_alleles.txt ./
    """
}




// ✅ **Generate Peptides**
process GeneratePeptides {
    tag "GeneratePeptides"

    output:
    path "mhc1_rand_peptides.fasta", emit: mhc1_fasta
    path "mhc2_rand_peptides.fasta", emit: mhc2_fasta

    script:
    """
    Rscript -e '
    library(stringi)

    generate_random_peptides <- function(num_peptides, peptide_length) {
      amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                       "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
      peptides <- replicate(num_peptides, stri_flatten(sample(amino_acids, peptide_length, replace = TRUE)))
      return(peptides)
    }

    set.seed(123)
    mhc1_peptides <- generate_random_peptides(1000000, 9)
    mhc2_peptides <- generate_random_peptides(1000000, 15)

    # ✅ Add FASTA headers to peptides
    mhc1_fasta <- paste0(">peptide_", seq_along(mhc1_peptides), "\\n", mhc1_peptides, collapse = "\\n")
    mhc2_fasta <- paste0(">peptide_", seq_along(mhc2_peptides), "\\n", mhc2_peptides, collapse = "\\n")

    writeLines(mhc1_fasta, "mhc1_rand_peptides.fasta")
    writeLines(mhc2_fasta, "mhc2_rand_peptides.fasta")
    '
    echo "✅ SUCCESS: Peptides generated and saved with correct FASTA headers."
    """
}

// ✅ **Split Peptides into Chunks**
process SplitPeptides {
    tag "SplitPeptides"

    input:
    path mhc1_fasta
    path mhc2_fasta

    output:
    path "mhc1_chunks/", emit: mhc1_chunks
    path "mhc2_chunks/", emit: mhc2_chunks

    script:
    """
    mkdir -p mhc1_chunks mhc2_chunks

    split -l 20000 ${mhc1_fasta} mhc1_chunks/mhc1_part_
    split -l 20000 ${mhc2_fasta} mhc2_chunks/mhc2_part_

    echo "✅ SUCCESS: Peptides split into chunks."
    """
}

process FilterPrecomputedAllelesFromRDS_MHC1 {
    publishDir "results/combined_mhc1", mode: 'copy'
    input:
    path allele_file
    path rds_file

    output:
    path "mhc1_alleles_filtered.txt", emit: filtered_alleles

    script:
    """
    Rscript -e '
    alleles <- readLines("${allele_file}")
    rds_data <- readRDS("${rds_file}")
    existing <- colnames(rds_data)
    missing <- setdiff(alleles, existing)
    writeLines(missing, "mhc1_alleles_filtered.txt")
    '
    """
}

process FilterPrecomputedAllelesFromRDS_MHC2 {
    publishDir "results/combined_mhc2", mode: 'copy'
    input:
    path allele_file
    path rds_file

    output:
    path "mhc2_alleles_filtered.txt", emit: filtered_alleles

    script:
    """
    Rscript -e '
    alleles <- readLines("${allele_file}")
    rds_data <- readRDS("${rds_file}")
    existing <- colnames(rds_data)
    missing <- setdiff(alleles, existing)
    writeLines(missing, "mhc2_alleles_filtered.txt")
    '
    """
}


process RunNetMHCpan {
    tag "NetMHCpan - ${sample_id}"
    errorStrategy 'ignore'


    input:
    tuple val(sample_id), path(mhc1_alleles), path(mhc1_chunks)

    output:
    path "mhc1_*.xls", emit: mhc1_results

    script:
    """
    echo "Running NetMHCpan for sample: ${sample_id}"
    for chunk in ${mhc1_chunks}/*; do
      chunk_name=\$(basename "\$chunk")

      grep -v '^\\s*\$' "${mhc1_alleles}" | xargs -I{} -P 7 bash -c '
        allele="{}"
        chunk="\$0"
        chunk_name="\$1"
        sample_id="\$2"
        output="mhc1_\${allele}_\${chunk_name}.xls"

        if [ ! -s "\$output" ]; then
          tcsh ${params.tool_mhc1} -a "\${allele}" -f "\${chunk}" -inptype 0 -l 9 -BA -xls -xlsfile "\${output}" > /dev/null 2>> error.log
        fi
      ' "\$chunk" "\$chunk_name" "${sample_id}"
    done

    echo "Finished NetMHCpan for sample: ${sample_id}"
    """
}



process RunNetMHCIIpan {
    errorStrategy 'ignore'

    tag "NetMHCIIpan - ${sample_id}"

    input:
    tuple val(sample_id), path(mhc2_alleles), path(mhc2_chunks)

    output:
    path "mhc2_*.xls", emit: mhc2_results

    script:
    """
    echo "Running NetMHCIIpan for sample: ${sample_id}"
    for chunk in ${mhc2_chunks}/*; do
      chunk_name=\$(basename "\$chunk")

      grep -v '^\\s*\$' "${mhc2_alleles}" | xargs -I{} -P 7 bash -c '
        allele="{}"
        chunk="\$0"
        chunk_name="\$1"
        sample_id="\$2"
        output="mhc2_\${allele}_\${chunk_name}.xls"

        if [ ! -s "\$output" ]; then

          tcsh ${params.tool_mhc2} -a "\${allele}" -f "\${chunk}" -inptype 0 -length 15 -xls -xlsfile "\${output}" > /dev/null 2>> error.log
        fi
      ' "\$chunk" "\$chunk_name" "${sample_id}"
    done

    echo "Finished NetMHCIIpan for sample: ${sample_id}"
    """
}

process CheckMissingChunks_MHC1 {
    input:
    path xls_files

    output:
    path "missing_chunks_mhc1.csv", emit: missing_chunks_mhc1

    script:
    """
    Rscript ${params.check_missing_script} . missing_chunks_mhc1.csv 10000
    """
}


process CheckMissingChunks_MHC2 {
    input:
    path xls_files

    output:
    path "missing_chunks_mhc2.csv", emit: missing_chunks_mhc2

    script:
    """
    Rscript ${params.check_missing_script} . missing_chunks_mhc2.csv 10000
    """
}

process RunNetMHCpan_Rerun {
    tag "Retry NetMHCpan - ${allele}_${chunk_name}"

    input:
    tuple val(allele), path(chunk_path), val(chunk_name), path(missing_chunks_path)

    output:
    path "mhc1_${allele}_*.xls", emit: rerun_results


    script:
    """
    export allele="${allele}"
    export chunk_name="${chunk_name}"
    export chunk_path="${chunk_path}"
    export missing_chunks_path="${missing_chunks_path}"

    output_file="mhc1_\${allele}_\${chunk_name}.xls"
    chunk_short=\$(echo "\${chunk_name}" | sed 's/^mhc1_part_//')

    should_run=\$(Rscript -e '
    df <- read.csv(Sys.getenv("missing_chunks_path"), stringsAsFactors = FALSE)
    allele <- Sys.getenv("allele")
    chunk <- Sys.getenv("chunk_name")
    chunk <- sub("^mhc1_part_", "", chunk)
    match <- any(df\$allele == allele & df\$chunk == chunk)
    cat(if (match) "TRUE" else "FALSE")
    ')

    if [ "\$should_run" = "TRUE" ]; then
        echo "🔁 Rerunning: allele=\$allele, chunk=\$chunk_name"
        tcsh ${params.tool_mhc1} -a "\$allele" -f "\$chunk_path" -inptype 0 -l 9 -BA -xls -xlsfile "\${output_file}" 2> error.log
        status=\$?

        if grep -q "not in allele list" error.log; then
            echo "⚠️ Skipped unsupported allele: \$allele"
            rm -f "\${output_file}" error.log
            exit 0
        fi

        if [ \$status -eq 0 ]; then
            echo "✅ Retry successful: \$allele, \$chunk_name"
        else
            echo "❌ Retry failed: \$allele, \$chunk_name, exit code \$status"
            cat error.log
        fi
    else
        echo "✔️ Chunk \$chunk_name for allele \$allele not flagged as missing, skipping."
    fi
    """



}

process RunNetMHCIIpan_Rerun {
    tag "Retry NetMHCIIpan - ${allele}_${chunk_name}"

    input:
    tuple val(allele), path(chunk_path), val(chunk_name), path(missing_chunks_path)

    output:
    path "mhc2_${allele}_*.xls", emit: rerun_results

    script:
    """
    output_file="mhc2_${allele}_${chunk_name}.xls"
    chunk_short="\$(echo "${chunk_name}" | sed 's/^mhc2_part_//')"

    # Use R to check if this allele + chunk are flagged as missing
    should_run=\$(Rscript -e '
      df <- read.csv("${missing_chunks_path}", stringsAsFactors = FALSE)
      allele <- "${allele}"
      chunk <- "${chunk_name}"
      chunk <- sub("^mhc2_part_", "", chunk)
      match <- any(df\$allele == allele & df\$chunk == chunk)
      cat(if (match) "TRUE" else "FALSE")
    ')

    if [ "\$should_run" = "TRUE" ]; then
        echo "🔁 Rerunning: allele=${allele}, chunk=${chunk_name}"

        tcsh ${params.tool_mhc2} -a "${allele}" -f "${chunk_path}" -inptype 0 -length 15 -xls -xlsfile "\${output_file}" 2> error.log
        status=\$?

        if grep -q "not in allele list" error.log; then
            echo "⚠️ Skipped unsupported allele: \${allele}"
            rm -f "\${output_file}" error.log
            exit 0
        fi

        if [ \$status -eq 0 ]; then
            echo "✅ Retry successful: \${allele}, \${chunk_name}"
        else
            echo "❌ Retry failed: \${allele}, \${chunk_name}, exit code \$status"
            cat error.log
        fi
    else
        echo "✔️ Chunk \${chunk_name} for allele \${allele} not flagged as missing, skipping."
    fi
    """
}
process UpdateGlobalMHC1Affinity {
    tag "Update global MHC-I affinity"
    publishDir "data/combined_mhc1_affinity.rds", mode: 'copy'

    input:
    path global_rds
    path new_rds

    output:
    path "${global_rds[0].getName()}", emit: updated_mhc1_rds


    script:
    """
    Rscript -e '
    message("🔁 Updating MHC-I global RDS...")
    global <- readRDS("${global_rds}")
    new <- readRDS("${new_rds}")

    new_cols <- setdiff(colnames(new), colnames(global))
    updated <- cbind(global, new[, new_cols, drop=FALSE])

    tmp_file <- "temp_mhc2_affinity.rds"
    final_file <- "${global_rds.getName()}"

    saveRDS(updated, tmp_file)

    if (!file.rename(tmp_file, final_file)) {
        stop("❌ Failed to move temporary file to final destination")
    }
    '
    """

}
process UpdateGlobalMHC2Affinity {
    tag "Update global MHC-II affinity"
    publishDir "data/combined_mhc2_affinity.rds", mode: 'copy'

    input:
    path global_rds
    path new_rds

    output:
    path "${global_rds[0].getName()}", emit: updated_mhc2_rds

    script:
    """
    Rscript -e '
      message("🔁 Updating MHC-II global RDS...")
      global <- readRDS("${global_rds}")
      new <- readRDS("${new_rds}")

      new_cols <- setdiff(colnames(new), colnames(global))
      updated <- cbind(global, new[, new_cols, drop=FALSE])

      saveRDS(updated, "${global_rds.getName()}")
    '
    """
}


process ParseNetMHCToRDS {
    tag "Parse NetMHC xls to RDS"

    input:
    path xls_files  // <- this is now a *list of files*
    path r_script

    output:
    path "*.rds", optional: true, emit: affinity_rds_mhc1


    script:
    """
    echo "Found \$(ls *.xls | wc -l) .xls files"
    Rscript  ${params.parse_netmhc_to_rds} . "mhc1_affinity_matrix.rds"
    """
}



process ParseNetMHCToRDS_II {
    tag "Parse NetMHC xls to RDS - II"

    input:
    path xls_dir
    path r_script

    output:
    path "*.rds", optional: true, emit: affinity_rds_mhc2

    script:
    """
    Rscript ${params.parse_netmhc_to_rds} ${xls_dir} ${xls_dir.simpleName}_affinity_matrix.rds
    """
}

process PrepareMGBSGenotypes {
    tag "Prepare MGBS Genotype RDS"

    input:
    path genotype_jsons


    output:
    path "mgbs_ready_genotypes.csv", emit: mgbs_rds_input

    script:
    """
    Rscript ${params.prepare_mgbs_genotypes} ${mapping_file} mgbs_ready_genotypes.rds ${genotype_jsons.join(" ")}
    """
}


process CalculateAllMGBS {
    tag "Calculate MGBS scores"
    publishDir "/home/dilara/projects/perimeta/results/nextflow_pipeline/MGBS/results/mgbs", mode: 'copy'

    input:
        tuple path(genotype_csv), path(mhc1_affinity_rds), path(mhc2_affinity_rds)

    output:
        path "*.rds", emit: mgbs_rds

    script:
    """
    Rscript ${params.calculate_mgbs_script} \\
      ${genotype_csv} \\
      ${mhc1_affinity_rds} \\
      ${mhc2_affinity_rds} \\
      MGBS_scores.rds \\
      ${params.kd_threshold}
    """
}



workflow {
    def peptide_out = GeneratePeptides()
    def chunk_out = SplitPeptides(peptide_out.mhc1_fasta, peptide_out.mhc2_fasta)

    def fastq_inputs

    if (params.data_format_genotype == "bam") {
        println "🔄 Using BAM input format"
        def bam_inputs = Channel
            .fromPath("${params.data_path}")
            .map { bam -> tuple(bam.getBaseName(), bam) }

        fastq_inputs = bam_inputs | BamConversion

    } else {
        println "📅 Using FASTQ input format"
        fastq_inputs = Channel
            .fromPath(params.data_path, type: 'dir')
            .filter { it.exists() }
            .map { dir ->
                def sample = dir.name
                def read1 = file("${dir}/${sample}_1.fastq.gz")
                def read2 = file("${dir}/${sample}_2.fastq.gz")
                if (!read1.exists() || !read2.exists()) return null
                tuple(sample, read1, read2)
            }
            .filter { it != null }
    }

    def genotype_jsons, parsed_csv_ch, mgbs_ready

    if (params.genotype_csv && file(params.genotype_csv).exists()) {
        println "📄 Using provided genotype CSV: ${params.genotype_csv}"
        parsed_csv_ch = Channel.value(file(params.genotype_csv))
        mgbs_ready = parsed_csv_ch
    } else {
        def hla_out
        if (params.hla_caller == "arcasHLA") {
            println "🧆 Running ArcasHLA"
            hla_out = fastq_inputs | ArcasHLA | PMappingAndFilter
            genotype_jsons = hla_out.map { it[1] }
            parsed_csv_ch = ParseGenotypesToCSV(genotype_jsons)
            mgbs_ready = PrepareMGBSGenotypes(genotype_jsons)
        } else if (params.hla_caller == "HLAHD") {
            println "🧆 Running HLA-HD"
            hla_out = fastq_inputs | HLAHD
            genotype_jsons = hla_out.map { it[1] }
            parsed_csv_ch = ParseHLAHDToCSV(genotype_jsons)
            mgbs_ready = parsed_csv_ch
        } else {
            error "❌ Unknown HLA caller selected: ${params.hla_caller}"
        }
    }

    def (mhc1_alleles, mhc2_alleles) = ExtractUniqueAlleles(parsed_csv_ch, file(params.convert_hla))

    if (params.use_precomputed_affinities) {
        println "🧠 Using precomputed affinities"

        def mhc1_rds_ch = Channel.value(file("/home/arne/projects/genomics_england/data/mhc1_rand_matrix.rds"))
        def mhc2_rds_ch = Channel.value(file("/home/arne/projects/genomics_england/data/mhc1_rand_matrix.rds"))

        def mhc1_filtered = FilterPrecomputedAllelesFromRDS_MHC1(mhc1_alleles, mhc1_rds_ch)
        def mhc2_filtered = FilterPrecomputedAllelesFromRDS_MHC2(mhc2_alleles, mhc2_rds_ch)

        def mhc1_combined = mhc1_filtered.map { a -> tuple("global", a) }.combine(chunk_out.mhc1_chunks)
        def mhc2_combined = mhc2_filtered.map { a -> tuple("global", a) }.combine(chunk_out.mhc2_chunks)

        def mhc1_results = mhc1_combined | RunNetMHCpan
        def mhc2_results = mhc2_combined | RunNetMHCIIpan

        def valid_mhc1_results = mhc1_results.filter { it.name.endsWith('.xls') && it.size() > 10 * 1024 }
        def valid_mhc2_results = mhc2_results.filter { it.name.endsWith('.xls') && it.size() > 10 * 1024 }

        def check_missing1 = CheckMissingChunks_MHC1(valid_mhc1_results.collect())
        def check_missing2 = CheckMissingChunks_MHC2(valid_mhc2_results.collect())

        def rerun_inputs_mhc1 = check_missing1.missing_chunks_mhc1
            .splitCsv(header: true)
            .map { row ->
                def chunk_name = "mhc1_part_${row.chunk}"
                def chunk_file = file("results/reference/mhc1_chunks/${chunk_name}")
                tuple(row.allele, chunk_file, chunk_name)
            }
            .combine(check_missing1.missing_chunks_mhc1)

        def rerun_inputs_mhc2 = check_missing2.missing_chunks_mhc2
            .splitCsv(header: true)
            .map { row ->
                def chunk_name = "mhc2_part_${row.chunk}"
                def chunk_file = file("results/reference/mhc2_chunks/${chunk_name}")
                tuple(row.allele, chunk_file, chunk_name)
            }
            .combine(check_missing2.missing_chunks_mhc2)

        def mhc1_rerun = rerun_inputs_mhc1 | RunNetMHCpan_Rerun
        def mhc2_rerun = rerun_inputs_mhc2 | RunNetMHCIIpan_Rerun

        // 📆 Wait for all MHC1 outputs (main + rerun), deduplicate by largest size
        def collected_primary = mhc1_results.collect()
        def collected_rerun = mhc1_rerun.collect()
        def all_mhc1_files = collected_primary
            .mix(collected_rerun)
            .flatten()
            .map { f -> tuple(f.getName(), f) }
            .groupTuple()
            .map { name, files -> files.max { it.size() } }
            .collect()


        def collected_mhc2_main = mhc2_results.collect()
        def collected_mhc2_rerun = mhc2_rerun.collect()

        def all_mhc2_results = collected_mhc2_main
            .mix(collected_mhc2_rerun)
            .flatten()
            .map { f -> tuple(f.getName(), f) }
            .filter { name, file -> file.size() > 10 * 1024 }
            .groupTuple()
            .map { name, files -> files.max { it.size() } }
            .collect()


        def parsed_mhc1 = ParseNetMHCToRDS(all_mhc1_files, file(params.parse_netmhc_to_rds)).affinity_rds_mhc1
        def parsed_mhc2 = ParseNetMHCToRDS_II(all_mhc2_results, file(params.parse_netmhc_to_rds)).affinity_rds_mhc2

        def updated_mhc1 = UpdateGlobalMHC1Affinity(file(params.mhc1_affinity), parsed_mhc1)
        def updated_mhc2 = UpdateGlobalMHC2Affinity(file(params.mhc2_affinity), parsed_mhc2)

        def mgbs_inputs = mgbs_ready.combine(updated_mhc1).combine(updated_mhc2)
        mgbs_inputs | CalculateAllMGBS

    } else {
        println "🧪 Calculating NetMHC affinities from scratch"

        def mhc1_combined = mhc1_alleles.map { a -> tuple("global", a) }.combine(chunk_out.mhc1_chunks)
        def mhc2_combined = mhc2_alleles.map { a -> tuple("global", a) }.combine(chunk_out.mhc2_chunks)
        // Run NetMHCpan
        def mhc1_results = mhc1_combined | RunNetMHCpan
        def mhc2_results = mhc2_combined | RunNetMHCIIpan

        // ✅ Filter only non-empty .xls files from NetMHCpan
        def valid_mhc1_results = mhc1_results.filter { it.name.endsWith('.xls') && it.size() > 10 * 1024 }

        // ✅ Only run CheckMissingChunks if we actually got some results
        def check_missing1 = valid_mhc1_results.ifEmpty { println "⚠️ No NetMHCpan results found — skipping missing chunk check"; return Channel.empty() }
                                            .collect()
                                            .ifNotEmpty { collected -> CheckMissingChunks_MHC1(collected) }

        // ✅ Always run CheckMissingChunks_MHC2 (NetMHCIIpan runs regardless)
        def check_missing2 = mhc2_results
            .ifEmpty { println "⚠️ No NetMHCIIpan results found — skipping missing chunk check"; return Channel.empty() }
            .collect()
            .ifNotEmpty { collected -> CheckMissingChunks_MHC2(collected) }



        def missing_csv1_path = check_missing1.missing_chunks_mhc1
        def missing_csv2_path = check_missing2.missing_chunks_mhc2

        def missing_csv_file = check_missing1.missing_chunks_mhc1
        def missing_csv_path = missing_csv_file.map { it.toAbsolutePath().toString() }

        def rerun_inputs_mhc1 = missing_csv_file
            .splitCsv(header: true)
            .combine(missing_csv_path)
            .map { row, path_str ->
                def chunk_name = "mhc1_part_${row.chunk}"
                def chunk_file = file("results/reference/mhc1_chunks/${chunk_name}")
                tuple(row.allele, chunk_file, chunk_name, path_str)
            }


        def rerun_inputs_mhc2 = missing_csv2_path
            .splitCsv(header: true)
            .map { row ->
                def chunk_name = "mhc2_part_${row.chunk}"
                def chunk_file = file("/results/reference/mhc2_chunks/${chunk_name}")
                tuple(row.allele, chunk_file, chunk_name, missing_csv2_path)
            }

        def mhc1_rerun = rerun_inputs_mhc1 | RunNetMHCpan_Rerun
        def mhc2_rerun = rerun_inputs_mhc2 | RunNetMHCIIpan_Rerun

        def all_mhc1_results = mhc1_results.mix(mhc1_rerun).collect()
        def all_mhc2_results = mhc2_results.mix(mhc2_rerun).collect()

        def mhc1_affinity = ParseNetMHCToRDS(all_mhc1_results, file(params.parse_netmhc_to_rds)).affinity_rds_mhc1
        def mhc2_affinity = ParseNetMHCToRDS_II(all_mhc2_results, file(params.parse_netmhc_to_rds)).affinity_rds_mhc2

        def mgbs_inputs = mgbs_ready.combine(mhc1_affinity).combine(mhc2_affinity)
        mgbs_inputs | CalculateAllMGBS
    }
}
