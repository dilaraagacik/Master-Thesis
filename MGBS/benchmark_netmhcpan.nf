#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow {
    def chunks = SplitPeptides(GeneratePeptides())

    chunks.mhc1_chunks.map { chunk -> tuple("benchmark", params.allele_mhc1, chunk) } | BenchmarkNetMHCpan
    chunks.mhc2_chunks.map { chunk -> tuple("benchmark", params.allele_mhc2, chunk) } | BenchmarkNetMHCIIpan
}

process GeneratePeptides {
    output:
    path "mhc1_rand_peptides.fasta", emit: mhc1_fasta
    path "mhc2_rand_peptides.fasta", emit: mhc2_fasta

    script:
    """
    Rscript -e '
    aa <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    writeLines(paste0(">peptide_", 1:1e6, "\n", replicate(1e6, paste(sample(aa, 9, TRUE), collapse=""))), "mhc1_rand_peptides.fasta")
    writeLines(paste0(">peptide_", 1:1e6, "\n", replicate(1e6, paste(sample(aa, 15, TRUE), collapse=""))), "mhc2_rand_peptides.fasta")
    '
    """
}

process SplitPeptides {
    input:
    path mhc1_fasta
    path mhc2_fasta

    output:
    path "mhc1_chunks/", emit: mhc1_chunks
    path "mhc2_chunks/", emit: mhc2_chunks

    script:
    """
    mkdir -p mhc1_chunks mhc2_chunks
    split -l 20000 ${mhc1_fasta} mhc1_chunks/part_
    split -l 20000 ${mhc2_fasta} mhc2_chunks/part_
    """
}

process BenchmarkNetMHCpan {
    tag "MHC-I"

    input:
    tuple val(sample_id), val(allele), path(chunk_dir)

    output:
    path("mhc1_all_chunks_${allele}.txt")

    script:
    """
    allele=${allele}
    for chunk in ${chunk_dir}/*; do
        chunk_name=\$(basename "\$chunk")
        tcsh ${params.tool_mhc1} -a ${allele} -f "\$chunk" -inptype 0 -l 9 -BA -xls -xlsfile mhc1_${allele}_\${chunk_name}.xls
    done
    echo "✅ NetMHCpan completed for allele ${allele}" > mhc1_all_chunks_\${allele}.txt
    """
}

process BenchmarkNetMHCIIpan {
    tag "MHC-II"

    input:
    tuple val(sample_id), val(allele), path(chunk_dir)

    output:
    path("mhc2_all_chunks_${allele}.txt")

    script:
    """
    allele=${allele}
    for chunk in ${chunk_dir}/*; do
        chunk_name=\$(basename "\$chunk")
        tcsh ${params.tool_mhc2} -a ${allele} -f "\$chunk" -inptype 0 -length 15 -xls -xlsfile mhc2_${allele}_\${chunk_name}.xls
    done
    echo "✅ NetMHCIIpan completed for allele ${allele}" > mhc2_all_chunks_${allele}.txt
    """
}