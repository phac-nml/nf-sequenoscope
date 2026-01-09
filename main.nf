#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { FILTER_ONT } from './modules/filter_ONT'
include { ANALYZE }    from './modules/analyze'
include { PLOT }       from './modules/plot'



workflow {
    if (params.help) { helpMessage(); return }

    def CLI_command = args[0]
    def run_batch = params.input_batch_tsv ? true : false
    ch_inputs = (params.input_batch_tsv) ? get_batch_inputs(params.input_batch_tsv) : get_single_inputs()

     // If batch, continue to Analyze; if single, stop here.
    if (run_batch) {
            ch_inputs = get_batch_inputs(params.input_batch_tsv)

            ch_inputs.branch {
                to_filter: it[2] != [] && !it[0].is_paired
                skip_filter: true
            }.set { ch_split }

            ch_split.to_filter.view { meta, fq, sum, ref, min, max ->
                "DEBUG [TO_FILTER]: ID=${meta.id} | Paired=${meta.is_paired} | Summary=${sum.name}"
            }

            ch_split.skip_filter.view { meta, fq, sum, ref, min, max ->
                "DEBUG [SKIP_FILTER]: ID=${meta.id} | Paired=${meta.is_paired} | Summary=${sum.name}"
            }
            
            ch_split.to_filter.dump(tag: 'TO_FILTER')
            ch_split.skip_filter.dump(tag: 'SKIP_FILTER')
            
            ch_filtered = FILTER_ONT(ch_split.to_filter.map { 
                meta, fastq, summary, ref, min, max -> 
                   return [ meta, fastq, summary, min, max ]
                }).results

            ch_filtered.view { "FILTER DEBUG: meta=${it[0]} | filtered_fq=${it[1]}" }
            ch_to_analyze_ont = ch_filtered
            .join(ch_inputs, by: 0) // Joins on the first element (meta)
            .map { 
                // meta           = index 0
                // filtered_fastq = index 1 (from FILTER_ONT)
                // orig_fastq     = index 2 (from ch_inputs original fastq)
                // summary        = index 3
                // reference      = index 4
                // min_ch         = index 5
                // max_ch         = index 6
                meta, filtered_fastq, orig_fastq, summary, reference, min_ch, max_ch ->

                println "----------------------------------------------------"
                println "ANALYZE PREP | ID: ${meta.id}"
                println "  - Modified File: ${filtered_fastq.name}"
                println "  - Metadata Ref:  ${reference.name ?: 'None'}"
                println "  - Channel Range: ${min_ch ?: 'All'}-${max_ch ?: 'All'}"
                println "----------------------------------------------------"

                return [ meta, filtered_fastq, summary, reference, min_ch, max_ch]
            }

            ch_to_analyze = ch_to_analyze_ont.mix(ch_split.skip_filter)

            ch_to_analyze.view { meta, fq, sum, ref, min, max ->
                "DEBUG: Sending to ANALYZE -> ID: ${meta.id}\n" +
                "  -> Fastq(s): ${fq}\n" +
                "  -> Summary: ${sum}\n" +
                "  -> Reference: ${ref}"
            }
            
            ch_analyzed = ANALYZE(ch_to_analyze).results
            ch_plot_input = ch_analyzed
                .map { meta, analyze_dir -> 
                    // Create a tuple where the barcode is the key
                    [ meta.barcode, meta, analyze_dir ] 
                }
                .groupTuple(by: 0) // Groups items by column barcode column together using index 0. That is meta.barcode
                .map { barcode, metas, dirs ->
                    // 2. Identify which directory is 'test' and which is 'control'
                    def test_idx = metas.findIndexOf { it.group == 'test' }
                    def ctrl_idx = metas.findIndexOf { it.group == 'control' }
                    
                    // Extract the correct directories
                    def test_dir = (test_idx != -1) ? dirs[test_idx] : dirs[0]
                    def ctrl_dir = (ctrl_idx != -1) ? dirs[ctrl_idx] : []
                    
                    // Create a new meta for the combined plot (e.g., "barcode1_comparison")
                    def plot_meta = [
                        id: "${barcode}_plots_comparison", 
                        barcode: barcode,
                        group: 'comparison'
                    ]
                    
                    return [ plot_meta, test_dir, ctrl_dir ]
                }
            // DEBUG VIEW: Print the paired channel content
            ch_plot_input.view { meta, t_dir, c_dir ->
                "DEBUG: Sending to PLOT -> Barcode: ${meta.barcode}\n" +
                "  -> Test Dir: ${t_dir}\n" +
                "  -> Ctrl Dir: ${c_dir ?: 'NONE (Single sample mode)'}"
            }    
            PLOT(ch_plot_input)
    }

    // SINGLE FILE MODE FILTERING MODE
    if (CLI_command == 'filter_ONT') {
        
        ch_inputs.view { meta, fastq, sum, ref, min, max ->
            "RUNNING: Starting FILTER_ONT for sample ${meta.id} ${sum} ${ref} ${min} ${max}"
        }

        //Need to map the 6-item ch_inputs to the 5-item tuple FILTER_ONT expects
        ch_filtered = FILTER_ONT(ch_inputs.map { meta, fq, sum, ref, min, max ->
            return [ meta, fq, sum, min, max ] }).results
                
        ch_filtered.view { meta, fastq -> 
            return """
            ------------------------------------------------------------------
            SUCCESS: FILTER_ONT for sample ${meta.id} is complete.
            - Output File: ${fastq.name}
            - Results moved to: ${params.output}/${meta.id}/filter_ONT/
            ------------------------------------------------------------------
            """.stripIndent()
        }
       
    }

    // SINGLE FILE MODE ANALYZE
    else if (CLI_command == 'analyze') {
        ch_inputs.view { meta, fastq, sum, ref, min, max ->
            "RUNNING: Starting ANALYZE for sample ${meta.id} ${fastq}"
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

    // SINGLE FILE MODE PLOT
    else if (CLI_command == 'plot') {
        // Handle input channel creation
        if (params.input_batch_tsv) {
            ch_plot_input = Channel.fromPath(params.input_batch_tsv)
                .splitCsv(sep: '\t', header: true)
                .map { row ->
                    def meta = [id: row.sample, group: row.group, barcode: row.barcode]
                    return [ meta, file(row.test_dir), row.control_dir ? file(row.control_dir) : [] ]
                }
        } 
        else if (params.test_dir && params.control_dir) {
            def meta = [id: params.output_prefix ?: "manual_plot"]
            ch_plot_input = Channel.of([
                meta, 
                file(params.test_dir), 
                params.control_dir ? file(params.control_dir) : [] 
            ])
        } 
        else {
            error "Plot mode requires either --input_batch_tsv or --test_dir and --control_dir"
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
      --force             Overwrite existing output files (default: ${params.force})

    Samplesheet Columns (CSV):
      sample              Unique ID for the run (used for file naming)
      fastq               Path to input FASTQ file
      summary             (Optional) Path to sequencing_summary.txt for filtering
      reference           Path to reference FASTA file
      min_ch, max_ch      (Optional) Pore channel range for subsetting
      group               (Optional) Specify two categories ('test' or 'control') for comparison plotting
      barcode             (Optional) Common ID to pair test/control samples for plotting

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



def get_batch_inputs(tsv) {
    def tsv_parent = file(tsv).parent // Get the parent directory of the TSV file
    return Channel.fromPath(tsv)
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def meta = [
                id: row.sample, 
                group: row.group ?: 'none', 
                barcode: row.barcode ?: row.sample,
                is_paired: row.fastq2 ? true : false
            ]
            def fastq_files = row.fastq2 ? [ tsv_parent.resolve(row.fastq.trim()), tsv_parent.resolve(row.fastq2.trim()) ] : [ tsv_parent.resolve(row.fastq.trim()) ]
            return [ meta, fastq_files, 
            row.sequence_summary_file ? tsv_parent.resolve(row.sequence_summary_file.trim()) : [], 
            row.reference_file ? tsv_parent.resolve(row.reference_file.trim()) : [], 
            row.min_ch ? row.min_ch.trim().toInteger() : null, 
            row.max_ch ? row.max_ch.trim().toInteger() : null ]
        }
}

def get_single_inputs() {
    if (!params.input_fastq) { error "Error: Must provide --input_fastq for this mode." }
    def meta = [id: params.output_prefix, group: 'none', barcode: params.output_prefix]
    def fastq_files = params.input_fastq2 ? 
        [ file(params.input_fastq), file(params.input_fastq2) ] : 
        file(params.input_fastq)
    println "DEBUG: Single input mode with FASTQ files: ${fastq_files}"
    return Channel.of([ meta, fastq_files, params.input_summary ? file(params.input_summary) : [], params.input_reference ? file(params.input_reference) : [], params.minimum_channel, params.maximum_channel ])
}
  
  
workflow.onComplete {
    def stats = workflow.stats
    def totalCpuHours = stats.succeedDuration.toHours() * stats.peakCpus 
    def msg = """
    ------------------------------------------------------------------
    PIPELINE SUMMARY
    ------------------------------------------------------------------
    Finished at  : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    Work Dir     : ${workflow.workDir}
    Exit status  : ${workflow.exitStatus}
    ------------------------------------------------------------------
    """.stripIndent()

    println msg
}