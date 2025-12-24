#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { FILTER_ONT } from './modules/filter_ONT'
include { ANALYZE }    from './modules/analyze'
include { PLOT }       from './modules/plot'

workflow {
    // BRANCH 1: FILTERING MODE
    if (params.mode == 'filter') {
        ch_inputs = (params.input_batch_tsv) ? get_batch_inputs(params.input_batch_tsv) : get_single_inputs()
        ch_inputs.view { meta, fastq, sum, ref, min, max ->
            "RUNNING: Starting FILTER_ONT for sample ${meta.id} [Mode: ${params.mode}]"
        }

        ch_filtered = FILTER_ONT(ch_inputs).results
        
        ch_filtered.view { meta, fastq, ref -> 
            return "SUCCESS: FILTER_ONT for sample ${meta.id} is complete. Results moved to: ${params.output}/${meta.id}/filter_ONT/"
        }
        
        // If batch, continue to Analyze; if single, stop here.
        if (run_batch) {
            ch_analyzed = ANALYZE(ch_filtered).results
            prepare_plot(ch_analyzed)
        }
    }

    // BRANCH 2: ANALYZE
    else if (params.mode == 'analyze') {
        ch_inputs.view { meta, fastq, sum, ref, min, max ->
            "RUNNING: Starting ANALYZE for sample ${meta.id}"
        }

        ch_analyzed = ANALYZE(ch_inputs).results
        ch_analyzed.view { meta, results_files ->
            "SUCCESS: ANALYZE for sample ${meta.id} is complete. Results moved to: ${params.output}/${meta.id}/analyze/"
        }
        
        // If batch, continue to Plot; if single, stop here.
        if (run_batch) {
            prepare_plot(ch_analyzed)
        }
    }

    // BRANCH 3: PLOT
    else if (params.mode == 'plot') {
        // Handle input channel creation
        if (params.input_batch_tsv) {
            ch_plot_input = Channel.fromPath(params.input_batch_tsv)
                .splitCsv(sep: '\t', header: true)
                .map { row ->
                    def meta = [id: row.sample, group: row.group, barcode: row.barcode]
                    return [ meta, file(row.test_dir), row.control_dir ? file(row.control_dir) : [] ]
                }
        } 
        else if (params.test_dir) {
            def meta = [id: params.output_prefix ?: "manual_plot"]
            ch_plot_input = Channel.of([
                meta, 
                file(params.test_dir), 
                params.control_dir ? file(params.control_dir) : [] 
            ])
        } 
        else {
            error "Plot mode requires either --input_batch_tsv or --test_dir"
        }

        // 1. Print that you are starting the plot mode
        ch_plot_input.view { meta, test, ctrl ->
            "RUNNING: Starting PLOT mode for sample ${meta.id}"
        }

        // 2. ACTUALLY CALL THE MODULE (This was missing)
        PLOT(ch_plot_input)

        // 3. Optional: View the success message
        PLOT.out.results.view { meta, plots ->
            "SUCCESS: PLOT complete for ${meta.id}. Results: ${params.output}/${meta.id}/plots/"
        }
    }

}

def helpMessage() {
    log.info"""
    nf-sequenoscope: Scalable Nextflow wrapper for Sequenoscope
    ==========================================================
    Usage:
      nextflow run main.nf --input samplesheet.csv [options]

    Main Options:
      --input             Path to CSV samplesheet (required for batch processing)
      --output            Directory to save results (default: ${params.output})
      --mode              Direct module override: [filter|analyze|plot]
      --force             Overwrite existing output files (default: ${params.force})

    Samplesheet Columns (CSV):
      sample              Unique ID for the run (used for file naming)
      fastq               Path to input FASTQ file
      summary             (Optional) Path to sequencing_summary.txt for filtering
      reference           Path to reference FASTA file
      min_ch, max_ch      (Optional) Pore channel range for subsetting
      group               (Optional) 'test' or 'control' for comparison plotting
      barcode             (Optional) Common ID to pair test/control runs

    Logic Overview:
      - If 'summary' and 'min_ch' are provided, the pipeline runs: 
        Filter_ONT -> Analyze.
      - If 'summary' is missing, the pipeline runs: 
        Analyze only.
      - If 'group' (test/control) and 'barcode' are provided, the pipeline:
        Pairs the analyzed results and runs the Plot module.

    Performance:
      -profile conda      Use automated Conda environment creation
      --threads           Number of threads per analysis job (default: ${params.threads})
    """.stripIndent()
}

// Show help message if --help is used
if (params.help) {
    helpMessage()
    exit 0
}



def get_batch_inputs(tsv) {
    return Channel.fromPath(tsv)
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def meta = [id: row.sample, group: row.group ?: 'none', barcode: row.barcode ?: row.sample]
            return [ meta, file(row.fastq), row.summary ? file(row.summary) : [], row.reference ? file(row.reference) : [], row.min_ch ?: 1, row.max_ch ?: 512 ]
        }
}

def get_single_inputs() {
    if (!params.input_fastq) { error "Error: Must provide --input_fastq for this mode." }
    def meta = [id: params.output_prefix, group: 'none', barcode: params.output_prefix]
    return Channel.of([ meta, file(params.input_fastq), params.input_summary ? file(params.input_summary) : [], params.input_reference ? file(params.input_reference) : [], params.minimum_channel, params.maximum_channel ])
}
  
  
