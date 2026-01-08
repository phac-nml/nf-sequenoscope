# nf-sequenoscope

This Nextflow wrapper extends Sequenoscope for high-throughput environments, enabling the parallel processing of multiple sequencing runs at scale. While adaptive sampling creates unique challenges in data management, this pipeline automates the orchestration across three specialized modules, ensuring consistent and reproducible interpretation across large-scale datasets. This pipeline addresses these through three specialized modules, using automated logic to detect and trigger the correct execution path based on your specific parameter inputs:

## Analyze Module

- **Functionality**: Acts as the core component, taking input from FASTQ files, a reference FASTA file, and an optional sequencing summary.
- **Operations**: Performs filtering, mapping, and generates comprehensive sequencing statistics, aiding in the initial processing and quality assessment of the sequencing data.

## Plot Module

- **Functionality**: Provides advanced visual analytics capabilities, turning complex data into interpretable visual formats.
- **Importance**: This module is essential for comparing and visualizing data from adaptive sampling experiments, using input from the analyze module to generate plots that highlight key differences and trends in the data.

## Filter_ONT Module

- **Focus**: This module focuses on the precise filtering of raw reads based on custom criteria such as channel number, duration, and quality scores, among others.
- **Utility**: It’s tailored for researchers needing detailed control over data inclusion in analyses, particularly useful in experiments with complex sampling strategies.

Together, these modules form a pipeline that not only simplifies the data processing workflow but also enhances the ability of researchers to conduct detailed and meaningful analyses of ONT sequencing data. The Nextflow wrapper ensures the pipeline is scalable, reproducible, and easily integrable into larger bioinformatics operations, making it an invaluable resource for genomic research.




nextflow run ../main.nf \
  --mode filter \
  --input_fastq ./fastq/barcode_1.fastq \
  --input_summary ./fastq/sequencing_summary_FAW13613_a5c2ca27_5cde2a43.txt \
  --minimum_q_score 1 \
  --minimum_channel 1 \
  --maximum_channel  256 \
  --output_prefix "barcode1_test"



nextflow run ../main.nf \
  --mode analyze \
  --threads 12 \
  --input_fastq ./sequenoscope_results/barcode1_test/filter_ONT/barcode1_test_filtered_fastq_subset.fastq \
  --input_reference ./fastq/Zymo_cleaned_ref.fasta \
  --output_prefix "barcode1_test"


  nextflow run ../main.nf \
  --mode plot \
  --test_dir ./sequenoscope_results/barcode1_test/analyze/ \
  --control_dir  ./sequenoscope_results/barcode1_test/analyze/  \
  --output_prefix "barcode1_test"



zip -r sequenoscope.zip . -x ".git/*" -x "conda-recipe/*" -x "mock_data/*" -x "./old_version/*" -x "./test_dir/*" -x "sequenoscope.egg-info/*" -x "./nextflow/*"


//batch mode
nextflow run ../main.nf  --mode filter --input_batch_tsv ./samplesheet.tsv  --output_prefix barcode1 --threads 12 -resume


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
seqtk sample -s100  /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/illumina_mock_data/SRA_Downloads/mock_community/ERR2984773_1.fastq 5000 > ERR2984773_1_EVEN.fastq 
seqtk sample -s100  /home/CSCScience.ca/kbessono/WORK/Sequenoscope/source/nf-sequenoscope/test_dir/fastq/illumina_mock_data/SRA_Downloads/mock_community/ERR2984773_2.fastq 5000 > ERR2984773_2_EVEN.fastq 