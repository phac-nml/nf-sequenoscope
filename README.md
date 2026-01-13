# nf-sequenoscope

This scalable Nextflow wrapper for Sequenoscope for high-throughput environments, enables parallel processing of multiple sequencing runs at scale. It automates orchestration across three specialized modules: Filter_ONT, Analyze, and Plot ensuring consistent and reproducible results generation across large datasets.

While adaptive sampling creates unique challenges in data management, this pipeline automates the orchestration across three specialized modules, ensuring consistent and reproducible interpretation across large-scale datasets. This pipeline addresses these through three specialized modules, using automated logic to detect and trigger the correct execution path based on your specific parameter inputs:

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
We used the ZymoBIOMICS mock communities to validate mapping and abundance estimation.

- Illumina EVEN/LOG: Data sourced from the [Loman Lab Mock Community project](https://github.com/LomanLab/mockcommunity) (Accessions: ERR2935805 and ERR2984773).

  - EVEN community: Contains 8 bacteria and 2 yeasts at equal genomic concentrations.

  - LOG community: Features the same organisms but with abundances distributed on a logarithmic scale.

- ONT Data: Subsampled Nanopore reads located in tests/data/ONT/ used for testing Adaptive Sampling filtration.

- References: A cleaned Zymo reference FASTA containing all taxa of interest is located in `tests/references/`.



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

* **Step 1: Filtration (ONT data only)**: If `min_ch`, `max_ch`, and a `sequence_summary_file` are provided, the pipeline automatically triggers the **FILTER_ONT** module to subset reads based on pore channel data.
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


