# Immune Selection Signal Analysis Pipeline

<img src="nextflow.png" alt="Nextflow" width="200"/>

A Nextflow pipeline to compute **immune selection signals** from HLA binding affinity in large cancer NGS datasets.

> The pipeline was developed and validated on the TCGA dataset and runs on HPC clusters (tested with Slurm). You may need to adapt some paths and configuration files for your cluster environment.

## Pipeline steps

1. [HLA Genotyping: arcasHLA / HLA-HD](https://github.com/RabadanLab/arcasHLA)
2. [MHC Binding Prediction: NetMHCpan / NetMHCIIpan](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)
3. [Neoantigen Depletion Analysis (Rsim)](https://pubmed.ncbi.nlm.nih.gov/31243365/)
4. MHC Genotype Binding Score (MGBS) calculation
5. PHBR score computation (Patient Harmonic Best Rank)
6. Immune cell deconvolution using [ConsensusTME](https://github.com/andrewGhazi/ConsensusTME)

## Usage

```bash
nextflow run main.nf --input /path/to/input --output /path/to/output
