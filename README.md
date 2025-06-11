# Immune Selection Signal Analysis Pipeline

<img src="nextflow.png" alt="Nextflow" width="200"/>

A Nextflow pipeline to compute **immune selection signals** from HLA binding affinity in large cancer NGS datasets.

> The pipeline was developed and validated on the TCGA dataset and runs on HPC clusters (tested with Slurm). You may need to adapt some paths and configuration files for your cluster environment.

## Pipeline steps

1. [HLA Genotyping: arcasHLA (v0.2.0)](https://github.com/RabadanLab/arcasHLA) / [HLA-HD (v1.3.0)][(https://w3.genome.med.kyoto-u.ac.jp/HLA-HD/)]
2. [MHC Binding Prediction: NetMHCpan (v4.0)](https://services.healthtech.dtu.dk/services/NetMHCpan-4.0/) / [NetMHCIIpan (v3.2)](https://services.healthtech.dtu.dk/services/NetMHCIIpan-3.2/)
3. Neoantigen Depletion Analysis (Rsim)
4. MHC Genotype Binding Score (MGBS) calculation
5. PHBR score computation (Patient Harmonic Best Rank)

## Usage

Example usage for PHBR
```bash
cd PHBR_Rsim
nextflow run main.nf -profile slurm --mode phbr --input_format maf
```
Example usage for MGBS
```
cd MGBS
nextflow run main.nf -profile slurm --use_precomputed_affinities true --genotypes path/to/data
```

### Available options

| Option                         | Description |
|--------------------------------|-------------|
| `--input_format`                | Mutation input file format (`vcf` or `maf`). |
| `--data_format_genotype`        | Input format for HLA genotyping (`bam` or `fastq`). |
| `--hla_caller`                  | Genotyping tool to use (`arcasHLA` or `hlahd`). |
| `--use_precomputed_affinities`  | Use precomputed affinity matrices (`true` or `false`). |
| `--mode`                        | Scoring module to run (`phbr` or `rsim`). |
| `--genotypes`                   | Path to CSV file with predefined HLA genotypes (optional). |
| `--kd_threshold`                | Threshold for predicted binding affinity in nM (default `500`). |
| `--barcode_map_csv`             | For TCGA dataset, CSV file mapping uuids to sample IDs and cancer types. |
| `--merged_counts`               | Path to precomputed gene expression matrix (optional, used in `rsim` mode). |

