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
    echo "Running HLA-HD for sample: ${sample_id}"

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
        echo "No result file found for sample: ${sample_id}" >&2
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

// Run ArcasHLA
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
        echo "Error: arcasHLA failed for ${sample_id}"
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




// Generate Peptides
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

    # Add FASTA headers to peptides
    mhc1_fasta <- paste0(">peptide_", seq_along(mhc1_peptides), "\\n", mhc1_peptides, collapse = "\\n")
    mhc2_fasta <- paste0(">peptide_", seq_along(mhc2_peptides), "\\n", mhc2_peptides, collapse = "\\n")

    writeLines(mhc1_fasta, "mhc1_rand_peptides.fasta")
    writeLines(mhc2_fasta, "mhc2_rand_peptides.fasta")
    '
    echo " SUCCESS: Peptides generated and saved with correct FASTA headers."
    """
}

// Split Peptides into Chunks
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

    echo "SUCCESS: Peptides split into chunks."
    """
}

process FilterPrecomputedAllelesFromRDS_MHC1 {

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

      if (length(missing) == 0) {
        writeLines("HLA-DONE", "mhc1_alleles_filtered.txt")
      } else {
        writeLines(missing, "mhc1_alleles_filtered.txt")
      }
    '
    """
}
process FilterPrecomputedAllelesFromRDS_MHC2 {

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

      if (length(missing) == 0) {
        writeLines("HLA-DONE", "mhc2_alleles_filtered.txt")
      } else {
        writeLines(missing, "mhc2_alleles_filtered.txt")
      }
    '
    """
}




process NetMHCpan {
    tag "NetMHCpan - ${allele}"
    publishDir "results/mhc1_xls", pattern: "mhc1_*.xls", mode: 'move'


    input:
    tuple val(allele), path(chunk_dir)

    output:
    path "mhc1_${allele}.xls", optional: true, emit: mhc1_result

    script:
    """
    set +e
    echo "Running NetMHCpan for allele: ${allele}"
    out_file="mhc1_${allele}.xls"
    found_any=0

    echo "📁 Using chunk directory: ${chunk_dir}"
    ls -lh "${chunk_dir}"

    for chunk in "${chunk_dir}"/mhc1_part_*; do
        if [ ! -f "\$chunk" ]; then
            echo "Skipping non-existent or non-regular file: \$chunk"
            continue
        fi

        echo "[\$(basename "\$chunk")] Running chunk for allele: ${allele}"
        tcsh ${params.tool_mhc1} -a "${allele}" -f "\$(realpath "\$chunk")" -inptype 0 -l 9 -BA -xls -xlsfile "temp.xls" > /dev/null 2> error.log
        status=\$?

        if grep -qi "could not find allele" error.log || grep -qi "not in allele list" error.log; then
            echo "Skipping unsupported allele: ${allele}"
            cat error.log
            rm -f temp.xls error.log
            continue
        fi

        if [[ \$status -eq 0 && -s temp.xls ]]; then
            echo "Successfully processed chunk: \$(basename "\$chunk")"
            cat temp.xls >> "\$out_file"
            found_any=1
        else
            echo "Chunk failed: \$(basename "\$chunk"), status=\$status"
            cat error.log
        fi

        rm -f temp.xls error.log
    done

    if [[ \$found_any -eq 0 ]]; then
        echo "No valid chunks processed for ${allele}"
        touch "mhc1_${allele}.xls"
    fi
    """
}


process NetMHCIIpan {
    tag "NetMHCIIpan - ${allele}"
    publishDir "results/mhc2_xls", pattern: "mhc2_*.xls", mode: 'copy'


    input:
    tuple val(allele), path(chunk_dir)

    output:
    path "mhc2_${allele}.xls", optional: true, emit: mhc2_result

    script:
    """
    set +e
    echo "Running NetMHCIIpan for allele: ${allele}"
    out_file="mhc2_${allele}.xls"
    found_any=0

    echo "Using chunk directory: ${chunk_dir}"
    ls -lh "${chunk_dir}"

    for chunk in "${chunk_dir}"/mhc2_part_*; do
        if [ ! -f "\$chunk" ]; then
            echo "Skipping non-existent or non-regular file: \$chunk"
            continue
        fi

        echo "[\$(basename "\$chunk")] Running chunk for allele: ${allele}"
        tcsh ${params.tool_mhc2} -a "${allele}" -f "\$(realpath "\$chunk")" -inptype 0 -length 15 -xls -xlsfile "temp.xls" > /dev/null 2> error.log
        status=\$?

        if grep -qi "could not find allele" error.log || grep -qi "not in allele list" error.log; then
            echo "Skipping unsupported allele: ${allele}"
            cat error.log
            rm -f temp.xls error.log
            # Don't exit! Just skip this chunk and continue
            continue
        fi


        if [[ \$status -eq 0 && -s temp.xls ]]; then
            echo "Successfully processed chunk: \$(basename "\$chunk")"
            cat temp.xls >> "\$out_file"
            found_any=1
        else
            echo "Chunk failed: \$(basename "\$chunk"), status=\$status"
            cat error.log
        fi

        rm -f temp.xls error.log
    done

    if [[ \$found_any -eq 0 ]]; then
        echo "No valid chunks processed for ${allele} (all skipped or unsupported)"
        touch "mhc2_${allele}.xls"
    fi

    """
}

process UpdateGlobalMHC1Affinity {
    tag "Update global MHC-I affinity"
    publishDir "data", mode: 'copy'

    input:
    path global_rds
    path new_rds

    output:
    path "mhc1_affinity.rds", emit: updated_mhc1_rds

    script:
    """
    Rscript -e '
      message("Updating MHC-I global RDS...")
      global <- readRDS("${global_rds}")
      new <- readRDS("${new_rds}")

      new_cols <- setdiff(colnames(new), colnames(global))

      if (length(new_cols) == 0 || nrow(new) == 0) {
        message("No new alleles or data rows to add — using existing RDS")
        file.copy("${global_rds}", "mhc1_affinity.rds")
      } else {
        updated <- cbind(global, new[, new_cols, drop=FALSE])
        saveRDS(updated, "mhc1_affinity.rds")
      }
    '
    """
}

process UpdateGlobalMHC2Affinity {
    tag "Update global MHC-II affinity"
    publishDir "data", mode: 'copy'

    input:
    path global_rds
    path new_rds

    output:
    path "mhc2_affinity.rds", emit: updated_mhc2_rds

    script:
    """
    Rscript -e '
      message("Updating MHC-II global RDS...")
      global <- readRDS("${global_rds}")
      new <- readRDS("${new_rds}")

      new_cols <- setdiff(colnames(new), colnames(global))

      if (length(new_cols) == 0 || nrow(new) == 0) {
        message("No new alleles or data rows to add — using existing RDS")
        file.copy("${global_rds}", "mhc2_affinity.rds")
      } else {
        updated <- cbind(global, new[, new_cols, drop=FALSE])
        saveRDS(updated, "mhc2_affinity.rds")
      }
    '
    """
}




process ParseNetMHCToRDS {
    tag "Parse NetMHC xls to RDS"

    input:
    path xls_dir
    path r_script

    output:
    path "*.rds", optional: true, emit: affinity_rds_mhc1

    script:
    """
    echo "Found \$(ls \${xls_dir}/*.xls | wc -l) .xls files"
    Rscript ${params.parse_netmhc_to_rds} ${xls_dir} mhc1_affinity_matrix.rds
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
    echo "Parsing MHC-II .xls files from directory: ${xls_dir}"
    ls -lh ${xls_dir}/*.xls || echo "No .xls files found"

    Rscript ${params.parse_netmhc_to_rds} ${xls_dir} mhc2_affinity_matrix.rds
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

def readAlleles = { file_ch ->
    file_ch
        .flatMap { file ->
            file.readLines().findAll { it?.trim() }
        }
}


def createAlleleChunkJobs = { alleles_ch, chunks_dir_ch ->
    chunks_dir_ch
        .map { dir -> file(dir).listFiles().findAll { it.name.startsWith('mhc') && it.size() > 0 } }
        .flatten()
        .combine(alleles_ch)
        .map { chunk, allele -> tuple(allele, chunk) }
}



workflow {
    def peptide_out = GeneratePeptides()
    def chunk_out = SplitPeptides(peptide_out.mhc1_fasta, peptide_out.mhc2_fasta)

    def fastq_inputs

    if (params.data_format_genotype == "bam") {
        println "Using BAM input format"
        def bam_inputs = Channel
            .fromPath("${params.data_path}")
            .map { bam -> tuple(bam.getBaseName(), bam) }

        fastq_inputs = bam_inputs | BamConversion

    } else {
        println "Using FASTQ input format"
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

    if (params.genotypes && file(params.genotypes).exists()) {
        println "Using provided genotype CSV: ${params.genotypes}"
        parsed_csv_ch = Channel.value(file(params.genotypes))
        mgbs_ready = parsed_csv_ch
    } else {
        def hla_out
        if (params.hla_caller == "arcasHLA") {
            println "Running ArcasHLA"
            hla_out = fastq_inputs | ArcasHLA | PMappingAndFilter
            genotype_jsons = hla_out.map { it[1] }
            parsed_csv_ch = ParseGenotypesToCSV(genotype_jsons)
            mgbs_ready = PrepareMGBSGenotypes(genotype_jsons)
        } else if (params.hla_caller == "HLAHD") {
            println "Running HLA-HD"
            hla_out = fastq_inputs | HLAHD
            genotype_jsons = hla_out.map { it[1] }
            parsed_csv_ch = ParseHLAHDToCSV(genotype_jsons)
            mgbs_ready = parsed_csv_ch
        } else {
            error "Unknown HLA caller selected: ${params.hla_caller}"
        }
    }

    def (mhc1_alleles, mhc2_alleles) = ExtractUniqueAlleles(parsed_csv_ch, file(params.convert_hla))

    if (params.use_precomputed_affinities) {
        println "Using precomputed affinities"

        // Load static precomputed RDS matrices
        def mhc1_rds_ch = Channel.value(file("/home/arne/projects/mhc_survival/data/mhc1_rand_matrix.rds"))
        def mhc2_rds_ch = Channel.value(file("/home/arne/projects/mhc_survival/data/mhc2_rand_matrix.rds"))

        // Filter out alleles already in the RDS matrices
        def mhc1_filtered = FilterPrecomputedAllelesFromRDS_MHC1(mhc1_alleles, mhc1_rds_ch)
        def mhc2_filtered = FilterPrecomputedAllelesFromRDS_MHC2(mhc2_alleles, mhc2_rds_ch)

        def mhc1_filtered_alleles = readAlleles(mhc1_filtered)
        def mhc2_filtered_alleles = readAlleles(mhc2_filtered)

        // Get chunk directories (so we loop internally inside the process)
        def mhc1_chunk_dir = chunk_out.mhc1_chunks.distinct()
        def mhc2_chunk_dir = chunk_out.mhc2_chunks.distinct()



        // One job per allele, directory passed as input
        def mhc1_jobs = mhc1_filtered_alleles.combine(mhc1_chunk_dir)
        def mhc2_jobs = mhc2_filtered_alleles.combine(mhc2_chunk_dir)

        def mhc1_results = mhc1_jobs | NetMHCpan
        def mhc2_results = mhc2_jobs | NetMHCIIpan

        def mhc1_result_dir = mhc1_results.collect().map { file("results/mhc1_xls") }
        def mhc2_result_dir = mhc2_results.collect().map { file("results/mhc2_xls") }

        def parsed_mhc1 = ParseNetMHCToRDS(mhc1_result_dir, file(params.parse_netmhc_to_rds)).affinity_rds_mhc1
        def parsed_mhc2 = ParseNetMHCToRDS_II(mhc2_result_dir, file(params.parse_netmhc_to_rds)).affinity_rds_mhc2

        def updated_mhc1 = UpdateGlobalMHC1Affinity(file(params.mhc1_affinity), parsed_mhc1)
        def updated_mhc2 = UpdateGlobalMHC2Affinity(file(params.mhc2_affinity), parsed_mhc2)

        def mgbs_inputs = mgbs_ready.combine(updated_mhc1).combine(updated_mhc2)
        mgbs_inputs | CalculateAllMGBS


    } else {
        println "Calculating NetMHC affinities from scratch"
        def mhc1_dir = chunk_out.mhc1_chunks.map { it.getParent() }.distinct()
        def mhc2_dir = chunk_out.mhc2_chunks.map { it.getParent() }.distinct()

        // One job per allele, receiving the full chunk directory
        def mhc1_combined = mhc1_alleles.combine(mhc1_dir)
        def mhc2_combined = mhc2_alleles.combine(mhc2_dir)

        // Run NetMHC tools
        def mhc1_results = mhc1_combined | NetMHCpan
        def mhc2_results = mhc2_combined | NetMHCIIpan

        def valid_mhc1_results = mhc1_results.filter { it.name.endsWith('.xls') }
        def valid_mhc2_results = mhc2_results.filter { it.name.endsWith('.xls') }

        def all_mhc1_results = valid_mhc1_results
        def all_mhc2_results = valid_mhc2_results


        def mhc1_affinity = ParseNetMHCToRDS(all_mhc1_results, file(params.parse_netmhc_to_rds)).affinity_rds_mhc1
        def mhc2_affinity = ParseNetMHCToRDS_II(all_mhc2_results, file(params.parse_netmhc_to_rds)).affinity_rds_mhc2

        def mgbs_inputs = mgbs_ready.combine(mhc1_affinity).combine(mhc2_affinity)
        mgbs_inputs | CalculateAllMGBS
    }
}
