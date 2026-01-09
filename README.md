# nf-sequenoscope

This scalable Nextflow wrapper for Sequenoscope for high-throughput environments, enabling parallel processing of multiple sequencing runs at scale. It automates orchestration across three specialized modules: Filter_ONT, Analyze, and Plot ensuring consistent and reproducible results generation across large datasets.

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

## Usage
### Batch Mode (Recommended)
The most powerful way to run nf-sequenoscope. This mode reads a Samplesheet (TSV) and automatically determines if it needs to filter, analyze, or plot the data


```
nextflow run main.nf --input_batch_tsv ./tests/samplesheets/samplesheet_ont.tsv --output batch_results_exp1 --threads 4
```

#### Processing logic
The pipeline uses an automated, conditional execution path based on the metadata provided in the samplesheet:

* **Step 1: Filtration**: If `min_ch`, `max_ch`, and a `sequence_summary_file` are provided, the pipeline automatically triggers the **FILTER_ONT** module to subset reads based on pore channel data.
* **Step 2: Core Analysis**: The **ANALYZE** module processes the reads (either filtered or raw) against the specified `reference_file` to generate mapping statistics and coverage profiles and other stats.
* **Step 3: Comparative Results Visualization**: Finally, the **PLOT** module identifies samples sharing a common `barcode` and pairs them based on the `group` field. This generates comparative plots for **Test** (e.g., Adaptive Sampling) vs. **Control** (e.g., Normal Run) conditions.


#### Samplesheet Format

| sample | fastq | fastq2 (Optional) | reference_file | group | barcode | 
| :--- | :--- | :--- | :--- | :--- | :--- | 
| sample1 | tests/data/reads1.fastq | | tests/data/ref.fa | test | barcode1 | 
| sample2 | tests/data/reads1.fastq | | tests/data/ref.fa | control | barcode1 |

>[!NOTE] 
>Paths in the samplesheet can be relative to the samplesheet file itself or use absolute paths

### Single-File Mode (Manual)
You can run individual modules for specific tasks for debugging or point testing purposes. 

#### Running Filter_ONT module only
Filters raw ONT reads based on channel number (e.g., to subset "adaptive sampled" vs "normal" reads from the same FASTQ input) and stops. Commands below will subset test and control samples based on the pore ranges.


```
nextflow run main.nf filter_ONT --input_fastq ./tests/data/ONT/barcode_1_small_2000.fastq --input_summary ./tests/data/ONT/test_sequencing_summary.txt --minimum_channel 1 --maximum_channel 256 --output results_filter_ONT --output_prefix test_sample

nextflow run main.nf filter_ONT --input_fastq ./tests/data/ONT/barcode_1_small_2000.fastq --input_summary ./tests/data/ONT/test_sequencing_summary.txt --minimum_channel 257 --maximum_channel 512 --output results_filter_ONT --output_prefix control_sample
```

#### Running ANALYZE module only
Runs the core Sequenoscope ANALYZE module on a single sample and stops. Here we would run ANALYZE on the subset reads from the FILTER_ONT module.


```
nextflow run main.nf analyze --input_fastq ./results_filter_ONT/control_sample/filter_ONT/control_sample_filtered_fastq_subset.fastq  --input_reference ./tests/references/Zymo_cleaned_ref.fasta  --output results_analyze --output_prefix test_sample --threads 4


nextflow run main.nf analyze --input_fastq ./results_filter_ONT/test_sample/filter_ONT/test_sample_filtered_fastq_subset.fastq  --input_reference ./tests/references/Zymo_cleaned_ref.fasta  --output results_analyze --output_prefix control_sample --threads 4

```

#### Running PLOT module only
Generates a comparison plots between two previously analyzed directories representing two conditions (e.g., adataptive sampling vs nomral sequencing).

```
nextflow run main.nf plot --test_dir ./results_analyze/test_sample/test_sample_analyze_results/  --control_dir ./results_analyze/control_sample/control_sample_analyze_results/  --output results_plot
```


nextflow run ../main.nf \
  filter_ONT \
  --input_fastq ./fastq/barcode_1.fastq \
  --input_summary ./fastq/sequencing_summary_FAW13613_a5c2ca27_5cde2a43.txt \
  --minimum_q_score 1 \
  --minimum_channel 1 \
  --maximum_channel  256 \
  --output_prefix "barcode1_test"



nextflow run ../main.nf \
   analyze \
  --threads 12 \
  --input_fastq ./sequenoscope_results/barcode1_test/filter_ONT/barcode1_test_filtered_fastq_subset.fastq \
  --input_reference ./fastq/Zymo_cleaned_ref.fasta \
  --output_prefix "barcode1_test"


  nextflow run ../main.nf \
  plot \
  --test_dir ./sequenoscope_results/barcode1_test/analyze/ \
  --control_dir  ./sequenoscope_results/barcode1_test/analyze/  \
  --output_prefix "barcode1_test"



zip -r sequenoscope.zip . -x ".git/*" -x "conda-recipe/*" -x "mock_data/*" -x "./old_version/*" -x "./test_dir/*" -x "sequenoscope.egg-info/*" -x "./nextflow/*" -x "./tests/*" -x ".*"


//batch mode
nextflow run ../main.nf  --input_batch_tsv ./samplesheet.tsv  --output_prefix barcode1 --threads 12 -resume


for i in {1..6}; do
    echo "Sampling 2000 reads from barcode_${i}..."
    seqtk sample -s100 "/home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/barcode_${i}.fastq" 2000 > "barcode_${i}_small_2000.fastq"
    sed -n '1~4p' barcode_${i}_small_2000.fastq | cut -d' ' -f1 | sed 's/^@//' >> sampled_ids.txt
done
echo "Total IDs collected: $(wc -l < sampled_ids.txt)"

head -n 1 /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/sequencing_summary_FAW13613_a5c2ca27_5cde2a43.txt > test_sequencing_summary.txt
grep -Ff sampled_ids.txt /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/sequencing_summary_FAW13613_a5c2ca27_5cde2a43.txt >> test_sequencing_summary.txt

echo "Total reads in new summary: $(($(wc -l < test_sequencing_summary.txt) - 1))"

#LOG
seqtk sample -s100  /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/illumina_mock_data/SRA_Downloads/mock_community/ERR2935805_1.fastq 50000 > ERR2935805_1_LOG.fastq 
seqtk sample -s100  /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/illumina_mock_data/SRA_Downloads/mock_community/ERR2935805_2.fastq 50000 > ERR2935805_2_LOG.fastq 

#EVEN
seqtk sample -s100  /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/illumina_mock_data/SRA_Downloads/mock_community/ERR2984773_1.fastq 50000 > ERR2984773_1_EVEN.fastq 
seqtk sample -s100  /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/illumina_mock_data/SRA_Downloads/mock_community/ERR2984773_2.fastq 50000 > ERR2984773_2_EVEN.fastq 


#single file mode

nextflow run ../main.nf \
  filter_ONT \
  --input_fastq /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/tests/data/ONT/barcode_1_small_2000.fastq \
  --input_summary ./fastq/sequencing_summary_FAW13613_a5c2ca27_5cde2a43.txt \
  --minimum_q_score 1 \
  --minimum_channel 1 \
  --maximum_channel  256 \
  --threads 12 \
  --output_prefix "barcode1_test"

  nextflow run ../main.nf \
  filter_ONT \
  --input_fastq /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/tests/data/ONT/barcode_1_small_2000.fastq \
  --input_summary ./fastq/sequencing_summary_FAW13613_a5c2ca27_5cde2a43.txt \
  --minimum_q_score 1 \
  --minimum_channel 257 \
  --maximum_channel  512 \
  --threads 12 \
  --output_prefix "barcode1_ctrl"

  nextflow run ../main.nf \
  analyze \
  --input_fastq /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/tests/data/Illumina/ERR2935805_1_LOG.fastq \
  --input_fastq2 /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/tests/data/Illumina/ERR2935805_2_LOG.fastq \
  --sequencing_type PE \
  --input_reference ./fastq/Zymo_cleaned_ref.fasta \
  --threads 12 \
  --output_prefix "IlluminaLOG"



  nextflow run ../main.nf  --input_batch_tsv ./samplesheets/samplesheet_illumina.tsv  --output_prefix illumina --threads 12 -resume

  nextflow run ../main.nf  --input_batch_tsv ./samplesheets/samplesheet_ont.tsv  --output_prefix ont --threads 12 -resume
