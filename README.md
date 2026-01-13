[![Powered by Sequenoscope](https://img.shields.io/badge/Powered%20by-Sequenoscope-blue)](https://github.com/phac-nml/sequenoscope)
![Nextflow](https://img.shields.io/badge/nextflow-%22%3E%3D23.04.0%22-brightgreen.svg)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/phac-nml/nf-sequenoscope)
![DSL2](https://img.shields.io/badge/Nextflow-DSL2-biueviolet)

## Contents <!-- omit from toc -->
- [About Sequenoscope](#about-sequenoscope)
- [Features](#features)
- [Module overview](#module-overview)
  - [Filter\_ONT Module](#filter_ont-module)
  - [Analyze module](#analyze-module)
  - [Plot Module](#plot-module)
- [Sample Data](#sample-data)
  - [ZymoBIOMICS Microbial Community Standards](#zymobiomics-microbial-community-standards)
- [Prerequisites \& Installation](#prerequisites--installation)
  - [Profiles](#profiles)
  - [Docker Permission Setup](#docker-permission-setup)
- [Usage](#usage)
  - [Batch Mode (Recommended for high throughput parallele processing)](#batch-mode-recommended-for-high-throughput-parallele-processing)
  - [Command Line Interface (CLI)](#command-line-interface-cli)
  - [Single-File Mode (Manual Walkthrough)](#single-file-mode-manual-walkthrough)
- [Citations](#citations)
- [Contacts](#contacts)

# nf-sequenoscope

This scalable [Nextflow](https://www.nextflow.io/) wrapper for **[Sequenoscope](https://github.com/phac-nml/sequenoscope)** designed for high-throughput environments, enables parallel processing of multiple sequencing runs at scale. It automates orchestration across three specialized modules: Filter_ONT, Analyze, and Plot ensuring consistent and reproducible results generation across large datasets.

While adaptive sampling creates unique challenges in data management, this pipeline automates the orchestration across three specialized modules, ensuring consistent and reproducible interpretation across large-scale datasets. This pipeline addresses these through three specialized modules, using automated logic to detect and trigger the correct execution path based on your specific parameter inputs:

## About Sequenoscope
This pipeline automates the orchestration of **[Sequenoscope](https://github.com/phac-nml/sequenoscope)**, a tool developed by the Public Health Agency of Canada (PHAC-NML). 

For detailed information on the underlying algorithms, mapping logic, and specific tool requirements, please refer to the original Sequenoscope repository at [https://github.com/phac-nml/sequenoscope](https://github.com/phac-nml/sequenoscope)

## Features
- End-to-End Automation in Batch Mode: Automatically chains filtering, analysis, and comparative plotting in a single run when using the batch mode.

- Batch Sample Processing: Handles complex manifests of multiple samples (ONT & Illumina) via a simple single TSV file.

- Flexible ONT Raw Reads Filtering: Specialized filtering for ONT data based on channel selection and other parameters allowing generation of sequencing data subsets from a single FASTQ file.

- Comparative Visualization (PLOT Module): Automatically identifies and pairs samples for Test vs. Control analysis—ideal for comparing Adaptive Sampling (AS) experiments against Normal/Standard sequencing runs to validate enrichment or depletion experimental conditions.

## Module overview
These modules form a pipeline that not only simplifies the data processing workflow but also enhances the ability of researchers to conduct detailed and meaningful analyses of ONT sequencing data. The Nextflow wrapper ensures the pipeline is scalable, reproducible, and easily integrable into larger bioinformatics operations, making it an invaluable resource for genomic research.


### Filter_ONT Module

- **Focus**: This module focuses on the precise filtering of raw reads based on custom criteria such as channel number, duration, and quality scores, among others.
- **Utility**: It’s tailored for researchers needing detailed control over data inclusion in analyses, particularly useful in experiments with complex sampling strategies.

### Analyze module

- **Functionality**: Acts as the core component, taking input from FASTQ files, a reference FASTA file, and an optional sequencing summary.
- **Operations**: Performs filtering, mapping, and generates comprehensive sequencing statistics, aiding in the initial processing and quality assessment of the sequencing data.

### Plot Module

- **Functionality**: Provides advanced visual analytics capabilities, turning complex data into interpretable visual formats.
- **Importance**: This module is essential for comparing and visualizing data from adaptive sampling experiments, using input from the analyze module to generate plots that highlight key differences and trends in the data.

## Sample Data
The repository includes test data to showcase and test pipeline performance across ONT and Illumina sequencing technologies and mock community compositions.

### ZymoBIOMICS Microbial Community Standards
We used the ZymoBIOMICS mock communities standards to validate mapping and abundance estimation.

- **Illumina EVEN/LOG Data:** Data sourced from the [Loman Lab Mock Community project](https://github.com/LomanLab/mockcommunity) (Accessions: ERR2935805 and ERR2984773).

  - **EVEN community:** Features 8 bacteria and 2 yeasts at equal genomic concentrations. The test set contains 50,000 paired-end reads (R1/R2).

  - **LOG community:** Features the same organisms with abundances distributed on a logarithmic scale. The test set contains 50,000 paired-end reads (R1/R2).

- **ONT Data:** Subsampled ONT reads from from LOG and EVEN ZymoBIOMICS mock communities standards for the 6 barcodes located in `tests/data/ONT/` folder. 
  - Each FASTQ file contains exactly 2,000 reads to facilitate rapid testing of the FILTER_ONT and ANALYZE modules.  
  -  A single `test_sequencing_summary.txt` (covering 12,000 reads) is provided, containing the shared 
  metadata for all 6 barcodes to enable testing of the FILTER_ONT module.
  - For the specific mapping of these 6 barcodes to their respective experimental conditions and channel ranges, please refer to the [manuscript](https://www.microbiologyresearch.org/content/journal/acmi/10.1099/acmi.0.001059.v2). 
  
**References:** A cleaned Zymo reference FASTA containing all taxa of interest is located in `tests/references/`.



## Prerequisites & Installation
This pipeline is built using Nextflow and requires it to be installed on your system.

- Nextflow: [Installation Guide](https://www.nextflow.io/docs/latest/install.html)
- Software Dependencies: You do not need to install Sequenoscope manually. The pipeline manages all dependencies through Conda, Docker, or Singularity. Choose the engine that best fits your environment.

### Profiles
Run the pipeline with the `-profile` to automatically handle all tool dependencies. If `-profile` is not specified the default `standard` profile is run.

| Profile | Description | Use Case Scenario |
| :--- | :--- | :--- |
| **standard** | Local execution | Software must be pre-installed in your `$PATH`. |
| **conda** | Conda-based environment | Personal workstations; no root access needed. |
| **docker** | Containerized (Docker) | Desktop/Linux servers; provides best reproducibility. |
| **singularity** | Containerized (Singularity) | Legacy HPC clusters; runs without root privileges. |
| **apptainer** | Containerized (Apptainer) | Modern HPC clusters; the successor to Singularity. |

### Docker Permission Setup
If you encounter a permission denied error when using `-profile docker` on Linux, your user likely lacks sudo access to the Docker daemon. Run these commands to fix it:

```
# Add your user to the docker group
sudo usermod -aG docker $USER
# Refresh your session to apply changes
newgrp docker
# Verify access (should show your local images)
docker images
```

## Usage
### Batch Mode (Recommended for high throughput parallele processing)
This is the most powerful way to run `nf-sequenoscope` pipeline. This mode reads a samplesheet (TSV) and automatically determines if it needs to filter or analyze the data first depending on available fields for each sample. All results for each run are saved in the `nf-sequenoscope-results` output folder by default (this could be customized by the `--output`).


```
# Run ONT batch (includes filtering if channel ranges are specified via min_ch and max_ch)
nextflow run main.nf --input_batch_tsv ./tests/samplesheets/samplesheet_ont.tsv --output batch_results_ont --threads 4

# Run Illumina batch (skips filtering automatically)
nextflow run main.nf --input_batch_tsv ./tests/samplesheets/samplesheet_illumina.tsv --output batch_results_illumina --threads 4
```

#### Processing logic
The pipeline uses an automated, conditional execution path based on the metadata provided in the samplesheet:

* **Step 1: Filtration (ONT data only)**: If `min_ch`, `max_ch`, and a `sequence_summary_file` fields are provided for a given sample, the pipeline automatically triggers the **FILTER_ONT** module to subset reads based on pore channel data.
* **Step 2: Core Analysis**: The **ANALYZE** module processes the reads (either filtered or raw) against the specified `reference_file` to generate mapping statistics and coverage profiles and other stats.
* **Step 3: Comparative Results Visualization**: Finally, the **PLOT** module identifies samples sharing a common `barcode` and pairs them based on the `group` field. This generates comparative plots for **Test** (e.g., Adaptive Sampling) vs. **Control** (e.g., Normal Run) conditions.
  
>[!NOTE]
>Illumina and ONT Data: For Illumina data or standard ONT runs (where adaptive sampling subsets are not required), the pipeline automatically skips **Step 1**. Since sequencing summaries and pore ranges do not apply to these data types, the execution begins directly with the ANALYZE module.

#### Samplesheet Format 
The samplesheet is a tab-separated (TSV) file used to drive the Batch Mode of the pipeline. It allows for the parallel processing of multiple runs, automatically determining whether to perform filtering or directly analyze samples ending with the comparative plotting based on the columns/fields provided.

Below is the example samplesheet with supported column names.

| sample | fastq | fastq2| reference_file | min_ch | max_ch | group | barcode | 
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | 
| sample1 | tests/data/reads1.fastq | | tests/data/ref.fa | 1 | 256 |test | barcode1 | 
| sample2 | tests/data/reads1.fastq | | tests/data/ref.fa | 257| 512 |control | barcode1 |

>[!NOTE] 
> - Paths in the samplesheet can be relative to the samplesheet file itself or use absolute paths
> - The `min_ch` and `max_ch` define the inclusive range of physical pore channels (e.g., 1–512 for MinION/GridION) to be extracted. These values must correspond exactly to the `channel` column in your `sequence_summary_file` to enable successful read sub-setting.
> - To ensure the PLOT module functions correctly, the `group` column strictly requires one of two specific values: `test` or `control`. Other values will not be recognized for sample pairing as this is essential for the PLOT module to perform comparative analysis.

#### Samplesheet Fields

The TSV samplesheet uses specific columns to orchestrate the workflow. Below is a description of each field and its impact on the pipeline logic

| **Column**                | **Description**                     | **Mandatory** | **Logic Impact**                                              |
| ------------------------- | ----------------------------------- | ------------- | ------------------------------------------------------------- |
| **sample**                | A unique identifier for the sample. | Yes           | Used for output directory naming and file prefixes.           |
| **fastq**                 | Path to the input FASTQ file.       | Yes           | Primary data source for all modules. If **fastq2** is not provided, assumes ONT long-read data.                         |
| **fastq2**                | Path to the R2 file for Illumina datasets.   | No            | If populated, the **ANALYZE** module runs in Paired-End mode |
| **reference_file**        | Path to the FASTA reference.        | Yes           | Used by **ANALYZE** for mapping and coverage stats.           |
| **sequence_summary_file** | ONT `sequencing_summary.txt`.       | No            | Required to run the **FILTER_ONT** module.                |
| **min_ch**                | Starting pore/channel (e.g. `1`).   | No            | Combined with `max_ch` to subset ONT reads for a given barcode.                   |
| **max_ch**                | Ending pore/channel (e.g. `256`).   | No            | Combined with `min_ch` to subset ONT reads for a given barcode.                   |
| **group**                 | Condition: `test` or `control`.     | No            | Required for the **PLOT** module to pair samples.             |
| **barcode**               | Common ID for pairing.              | No            | Samples with the same barcode are paired for comparison.      |


### Command Line Interface (CLI)
For a full list of all available parameters and their default values, run the pipeline with the `--help` flag.


```text
nf-sequenoscope: Scalable Nextflow wrapper for Sequenoscope
version: 1.0.0
==========================================================

1. BATCH MODE (Recommended)
    Automatically handles filtering, analysis, and plotting for multiple samples.
    Usage:
        nextflow run main.nf --input_batch_tsv samplesheet.tsv [options]

2. SINGLE MODULE MODE (Manual/Debugging)
    Run specific modules individually. Useful for troubleshooting.
    Usage:
        nextflow run main.nf <MODULE_NAME> [options]

   Available Modules:
     filter_ONT    - Subset reads based on pore channel ranges
     analyze       - Run alignment and mapping statistics
     plot          - Generate comparative enrichment plots

Main Options:
  --input_batch_tsv   Path to TSV samplesheet (required for batch processing)
  --output            Directory to save results (default: nf-sequenoscope-results)
  --threads           Number of threads per analysis job (default: 1)
  --force             Overwrite existing output files (default: false)

Module-Specific Options:

  [Filter_ONT]
  --input_fastq       Path to raw FASTQ
  --input_summary     ONT sequencing_summary.txt
  --minimum_channel   Min pore channel (default: 1)
  --maximum_channel   Max pore channel (default: 512)
  --minimum_q_score   Min Q-score threshold (default: 0)
  --minimum_length    Min read length (default: 0)

  [Analyze]
  --input_reference   Path to reference FASTA
  --sequencing_type   'SE' or 'PE' (default: SE)
  --quality_threshold Min map quality (default: 15)
  --minimap2_kmer     K-mer size for alignment (default: 15)

  [Plot]
  --test_dir          Path to Test results
  --control_dir       Path to Control results
  --time_bin_unit     Bin size for time plots (default: 15m)
  --adaptive_sampling Is this an AS run? (default: true)

Samplesheet Columns (TSV):
  sample              Unique ID for the run (used for file naming)
  fastq               Path to input FASTQ file
  fastq2              (Optional) Path to R2 for Illumina paired-end data
  sequence_summary_file (Optional) Path to ONT sequencing_summary.txt for filtering
  reference_file      Path to reference FASTA file
  min_ch, max_ch      (Optional) Pore channel range for subsetting (requires summary)
  group               Must be 'test' or 'control' for comparative plotting
  barcode             Common ID to pair 'test' and 'control' samples

Logic Overview:
  - Step 1 (Filter): Triggered if 'sequence_summary_file', 'min_ch', and 'max_ch' are present.
  - Step 2 (Analyze): Runs on all samples (filtered ONT, raw ONT, or Illumina).
  - Step 3 (Plot): Triggered for samples sharing a 'barcode' with 'test'/'control' groups.

Profiles:
  -profile conda      Managed Conda environments
  -profile docker     Docker containers (requires sudo usermod -aG docker $USER)
  -profile singularity/apptainer  HPC-friendly containers
```      

### Single-File Mode (Manual Walkthrough)
You can run individual modules for debugging or point testing purposes. Using the sample data stored in the tests/ folder, you can manually step through the entire pipeline—from filtering raw reads to generating final comparative plots.

>[!TIP] While manual mode is great for debugging, Batch Mode is highly recommended for large sample sizes as it automates these steps into a single command.

#### Step 1: Running Filter_ONT module (Filter raw ONT reads)
Filters raw ONT reads based on channel number (e.g., to subset "adaptive sampled" vs "normal" reads from the same FASTQ input) and stops. 


First, we subset "adaptive sampled" vs "normal" reads from the same FASTQ input based on physical pore ranges.


```
# Filter for 'Test' group undergoing AS sequencing (Channels 1-256)
nextflow run main.nf filter_ONT --input_fastq ./tests/data/ONT/barcode_1_small_2000.fastq --input_summary ./tests/data/ONT/test_sequencing_summary.txt --minimum_channel 1 --maximum_channel 256 --output results_filter_ONT --output_prefix test_sample

# Filter for 'Control' group undergoing Normal sequencing (Channels 257-512)
nextflow run main.nf filter_ONT --input_fastq ./tests/data/ONT/barcode_1_small_2000.fastq --input_summary ./tests/data/ONT/test_sequencing_summary.txt --minimum_channel 257 --maximum_channel 512 --output results_filter_ONT --output_prefix control_sample
```

#### Step 2: Running ANALYZE module (Analyze subsetted reads)
Runs the core Sequenoscope ANALYZE module on a single sample and stops. Here we would run ANALYZE on the subset reads from the FILTER_ONT module.

Next, run the core ANALYZE module on the outputs from Step 1 to generate read mapping statistics against the Zymo reference.

```
# Run Analyze on Test Sample
nextflow run main.nf analyze --input_fastq  ./results_filter_ONT/test_sample/filter_ONT/test_sample_test_sample_filtered_fastq_subset.fastq --input_reference ./tests/references/Zymo_cleaned_ref.fasta  --output results_analyze --output_prefix test_sample --threads 4

# Run Analyze on Control Sample
nextflow run main.nf analyze --input_fastq ./results_filter_ONT/control_sample/filter_ONT/control_sample_control_sample_filtered_fastq_subset.fastq  --input_reference ./tests/references/Zymo_cleaned_ref.fasta  --output results_analyze --output_prefix control_sample --threads 4

```

#### Step 3: Running PLOT module (Generate comparative plots)
Finally, use the PLOT module to compare the two experimental conditions (e.g., adataptive sampling vs nomral sequencing).

```
nextflow run main.nf plot --test_dir ./results_analyze/test_sample/test_sample_analyze_results/  --control_dir ./results_analyze/control_sample/control_sample_analyze_results/  --output results_plot
```


## Citations
If you use Sequenoscope in your research, please cite:

Meknas, A., Bessonov, K., Eagle, S. H., Peterson, C. L., Robertson, J., Ricker, N., ... & Reimer, A. (2025). **Sequenoscope: A Modular Tool for Nanopore Adaptive Sequencing Analytics and Beyond.** *Access Microbiology*, 001059-v2. https://doi.org/10.1099/acmi.0.001059.v2


## Contacts
- **Kyrylo Bessonov**: kyrylo.bessonov@phac-aspc.gc.ca 
- **Abdallah Meknas**: abdallahmeknas@gmail.com