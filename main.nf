#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Import modules
include { FILTER_ONT } from './modules/filter_ONT'
include { ANALYZE }    from './modules/analyze'
include { PLOT }       from './modules/plot'

workflow {
    def MODULE_FILTER_NAME = 'filter_ONT'
    def MODULE_ANALYZE_NAME = 'analyze'
    def MODULE_PLOT_NAME    = 'plot'

    if (params.help) { helpMessage(); exit 0}
    if (params.version) { println "nf-sequenoscope version: ${workflow.manifest.version}"; exit 0 }

    
    def CLI_command = args.size() > 0 ? args[0] : null

     // If batch, continue to Analyze; if single, stop here.
    if (params.input_batch_tsv) {
            ch_inputs = get_batch_inputs(params.input_batch_tsv)

            ch_inputs.branch {
                to_filter: it[2] != [] && !it[0].is_paired
                skip_filter: true
            }.set { ch_split }

            ch_split.to_filter.view{ meta, fq, sum, ref, min, max ->
                log.debug "DEBUG [TO_FILTER]: ID=${meta.id} | Paired=${meta.is_paired} | Summary=${sum.name}"
            }

            ch_split.skip_filter.view { meta, fq, sum, ref, min, max ->
                log.debug "DEBUG [SKIP_FILTER]: ID=${meta.id} | Paired=${meta.is_paired} | Summary=${sum.name}"
            }
            
            ch_split.to_filter.dump(tag: 'TO_FILTER')
            ch_split.skip_filter.dump(tag: 'SKIP_FILTER')
            
            ch_filtered = FILTER_ONT(ch_split.to_filter.map { 
                meta, fastq, summary, ref, min, max -> 
                   return [ meta, fastq, summary, min, max ]
                }).results

            ch_filtered.view { log.debug "FILTER DEBUG: meta=${it[0]} | filtered_fq=${it[1]}" }
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

                log.debug """
                ----------------------------------------------------
                ANALYZE CHANNEL PREP | ID: ${meta.id}
                - Modified File: ${filtered_fastq.name}
                - Metadata Ref:  ${reference.name ?: 'None'}
                - Channel Range: ${min_ch ?: 'All'}-${max_ch ?: 'All'}
                ----------------------------------------------------
                """.stripIndent()

                return [ meta, filtered_fastq, summary, reference, min_ch, max_ch]
            }

            ch_to_analyze = ch_to_analyze_ont.mix(ch_split.skip_filter)

            ch_to_analyze.view { meta, fq, sum, ref, min, max ->
                log.debug "DEBUG: Sending to ANALYZE -> ID: ${meta.id}\n" +
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
                log.debug "DEBUG: Sending to PLOT -> Barcode: ${meta.barcode}\n" +
                "  -> Test Dir: ${t_dir}\n" +
                "  -> Ctrl Dir: ${c_dir ?: 'NONE (Single sample mode)'}"
            }    
            PLOT(ch_plot_input)
    }

    // SINGLE FILE MODE FILTERING MODE
    else if (CLI_command == MODULE_FILTER_NAME) {
        ch_inputs = get_single_inputs() 
        ch_inputs.view { meta, fastq, sum, ref, min, max ->
            log.info "RUNNING: Starting FILTER_ONT for sample ${meta.id} ${sum} ${ref} ${min} ${max}"
        }

        //Need to map the 6-item ch_inputs to the 5-item tuple FILTER_ONT expects
        ch_filtered = FILTER_ONT(ch_inputs.map { meta, fq, sum, ref, min, max ->
            return [ meta, fq, sum, min, max ] }).results
                
        ch_filtered.view { meta, fastq -> 
            log.info """
            ------------------------------------------------------------------
            SUCCESS: FILTER_ONT for sample ${meta.id} is complete.
            - Output File: ${fastq.name}
            - Results moved to: ${params.output}/${meta.id}/filter_ONT/
            ------------------------------------------------------------------
            """.stripIndent()
        }
       
    }

    // SINGLE FILE MODE ANALYZE
    else if (CLI_command == MODULE_ANALYZE_NAME) {
        ch_inputs = get_single_inputs() 
        ch_inputs.view { meta, fastq, sum, ref, min, max ->
            log.info "RUNNING: Starting ANALYZE for sample ${meta.id} ${fastq}"
        }

        ch_analyzed = ANALYZE(ch_inputs).results
        ch_analyzed.view { meta, results_files ->
            log.info "SUCCESS: ANALYZE for sample ${meta.id} is complete. Results moved to: ${params.output}/${meta.id}/analyze/"
        }
        
    }

    // SINGLE FILE MODE PLOT
    else if (CLI_command == MODULE_PLOT_NAME) {
        
        if (params.test_dir && params.control_dir) {
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
            log.info"RUNNING: Starting PLOT mode for sample ${meta.id}"
        }

        // 2. ACTUALLY CALL THE MODULE (This was missing)
        PLOT(ch_plot_input)

        // 3. Optional: View the success message
        PLOT.out.results.view { meta, plots ->
            log.info"SUCCESS: PLOT complete for ${meta.id}. Results: ${params.output}/${meta.id}/plots/"
        }
    }
    else {
            log.error "ERROR: Invalid usage. Provide --input_batch_tsv OR a module name (filter_ONT|analyze|plot)."
            helpMessage()
            exit 1
    }

}





def helpMessage() {
    log.info"""
    nf-sequenoscope: Scalable Nextflow wrapper for Sequenoscope
    version: ${workflow.manifest.version}
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
      --output            Directory to save results (default: ${params.output})
      --threads           Number of threads per analysis job (default: ${params.threads})
      --force             Overwrite existing output files (default: ${params.force})

    Module-Specific Options:
    
      [Filter_ONT]
      --input_fastq       Path to raw FASTQ
      --input_summary     ONT sequencing_summary.txt
      --minimum_channel   Min pore channel (default: ${params.minimum_channel})
      --maximum_channel   Max pore channel (default: ${params.maximum_channel})
      --minimum_q_score   Min Q-score threshold (default: ${params.minimum_q_score})
      --minimum_length    Min read length (default: ${params.minimum_length})

      [Analyze]
      --input_reference   Path to reference FASTA
      --sequencing_type   'SE' or 'PE' (default: ${params.sequencing_type})
      --quality_threshold Min map quality (default: ${params.quality_threshold})
      --minimap2_kmer     K-mer size for alignment (default: ${params.minimap2_kmer})

      [Plot]
      --test_dir          Path to Test results
      --control_dir       Path to Control results
      --time_bin_unit     Bin size for time plots (default: ${params.time_bin_unit})
      --adaptive_sampling Is this an AS run? (default: ${params.adaptive_sampling})

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
      -profile docker     Docker containers (requires sudo usermod -aG docker \$USER)
      -profile singularity/apptainer  HPC-friendly containers
    """.stripIndent()
}



def get_batch_inputs(tsv) {
    def tsv_parent = file(tsv).parent // Get the parent directory of the TSV file
    def required_columns = [
        'sample', 'fastq', 'fastq2', 'sequence_summary_file', 
        'reference_file', 'min_ch', 'max_ch', 'group', 'barcode'
    ]
    return Channel.fromPath(tsv)
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            // 1. VALIDATE HEADERS (Check if any required column is missing from the Map keys)
            def missing_cols = required_columns.findAll { !row.containsKey(it) }
            if (missing_cols) {
                error """
                ------------------------------------------------------------------
                MISSING COLUMNS IN TSV
                The following required columns are missing: ${missing_cols.join(', ')}
                Check your TSV header names and ensure they are Tab-separated.
                ------------------------------------------------------------------
                """.stripIndent()
            }
            
            def fastq_files = row.fastq2 ? [ tsv_parent.resolve(row.fastq.trim()), tsv_parent.resolve(row.fastq2.trim()) ] : [ tsv_parent.resolve(row.fastq.trim()) ]
            def grp = row.group ? row.group.trim().toLowerCase() : ''
            def allowed_groups = ['test', 'control'] //only these two values are allowed
            // 2. VALIDATE 'group' column VALUES
            if ( !allowed_groups.contains(grp) ) {
                error """
                ------------------------------------------------------------------
                INVALID SAMPLESHEET ENTRY
                Sample: ${row.sample}
                Error : The 'group' column value '${row.group}' is invalid.
                Fix   : Use 'test', 'control'.
                ------------------------------------------------------------------
                """.stripIndent()
            }
            def meta = [
                id: row.sample, 
                group: row.group ?: '', 
                barcode: row.barcode ?: row.sample,
                is_paired: row.fastq2 ? true : false
            ]
            
            return [ meta, fastq_files, 
            row.sequence_summary_file ? tsv_parent.resolve(row.sequence_summary_file.trim()) : [], 
            row.reference_file ? tsv_parent.resolve(row.reference_file.trim()) : [], 
            row.min_ch ? row.min_ch.trim().toInteger() : null, 
            row.max_ch ? row.max_ch.trim().toInteger() : null ]
        }
}

def get_single_inputs() {
    if (!params.input_fastq) { error "Error: Must provide --input_fastq for this mode." }
    def meta = [id: params.output_prefix, group: '', barcode: params.output_prefix]
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