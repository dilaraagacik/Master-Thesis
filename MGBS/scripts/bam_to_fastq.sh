#!/usr/bin/env bash
set -euo pipefail

sample_id=$1
bam_file=$2
threads=$3

echo "🔃 Sorting and indexing full BAM..."
samtools sort -@ "$threads" -o "${sample_id}.sorted.bam" "$bam_file"
samtools index "${sample_id}.sorted.bam"

echo "🔄 Converting full BAM to FASTQ..."

samtools fastq -@ "$threads" \
    -1 "${sample_id}_1.fastq" \
    -2 "${sample_id}_2.fastq" \
    -0 /dev/null -s /dev/null -n \
    "${sample_id}.sorted.bam"

gzip "${sample_id}_1.fastq"
gzip "${sample_id}_2.fastq"

echo "✅ FASTQ files created from full BAM:"
echo "   → ${sample_id}_1.fastq.gz"
echo "   → ${sample_id}_2.fastq.gz"
