# A Nextflow pipeline to analyse HLA binding affinity based immune selection signals in cancer next generation sequencing data.

Immune Selection Signal Analysis Pipeline

This repository contains a Nextflow-based pipeline developed to analyze HLA binding affinity-based immune selection signals in large cancer NGS datasets. The pipeline automates the computation of three key immunogenomic metrics:

MHC Genotype Binding Score (MGBS)

Patient Harmonic Best Rank (PHBR)

Rsim metric for neoantigen depletion

These metrics enable the quantification of immune selection signals and their correlation with immune infiltration across different cancer types.

Features
Modular, reproducible pipeline implemented in Nextflow

Flexible input formats: supports RNA-seq and whole exome sequencing (WES) data

Automated workflows for:

HLA genotyping with arcasHLA and HLA-HD

MHC binding prediction with NetMHCpan and NetMHCIIpan

Neoantigen depletion analysis (Rsim)

Immune cell deconvolution using ConsensusTME

Supports large-scale pan-cancer analyses (validated on TCGA dataset)

Pipeline Overview
MGBS Workflow
Estimates the peptide-binding capacity of a patient’s HLA genotype using randomly generated peptides:

MGBS-I (MHC class I)

MGBS-II (MHC class II)

PHBR / Rsim Workflow
Predicts the immunogenic potential of somatic mutations and quantifies neoantigen depletion:

PHBR scores for MHC-I and MHC-II

Rsim metric to assess negative selection signals in HLA-binding regions

Installation
Requirements:

Nextflow

Conda (for environment management)

Slurm or other compatible workload manager (for HPC)

Clone the repository:

git clone https://github.com/dilaraagacik/Master-Thesis.git

Usage
Run MGBS pipeline:

nextflow run main.nf --input /path/to/input --output /path/to/output

Run PHBR + Rsim pipeline:


nextflow run main.nf --input /path/to/input --output /path/to/output
Options and configuration parameters can be set in nextflow.config.

Data Sources
TCGA WES + RNA-seq data

Reference HLA allele sequences (IPD-IMGT/HLA database)

NetMHCpan / NetMHCIIpan models

References
Claeys and Van den Eynden, 2024 - MGBS metric

Marty et al., 2017 - PHBR metric

Van den Eynden et al., 2019 - Rsim metric



