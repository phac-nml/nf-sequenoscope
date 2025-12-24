process RUN_ANALYZE {
    tag "${meta.id}"

    publishDir "${params.output}/${meta.id}/analyze", 
        mode: 'copy', 
        saveAs: { filename -> file(filename).name }

    input:
    tuple val(meta), path(fastqs), path(input_summary), path(reference), val(min_ch), val(max_ch)

    output:
    tuple val(meta), path("analyze_results/*"), emit: analysis_results

    script:
    // Simply list the parameter names that match the long-form flags
    def analyze_flags = [
        'minimum_coverage', 'threads', 'minimum_read_length', 
        'maximum_read_length', 'trim_front_bp', 'trim_tail_bp', 
        'quality_threshold', 'start_time', 'end_time', 
        'output_prefix', 'minimap2_kmer'
    ]

    def options = analyze_flags.collect { flag ->
        if (params[flag] != null && params[flag] != get_default_analyze(flag)) {
            return "--${flag} ${params[flag]}"
        }
    }.findAll().join(' ')

    if (params.force) options += " --force"
    
    // Handle the sequencing summary and type explicitly as they are often required
    def seq_sum_opt = (input_summary && !input_summary.empty()) ? "--sequencing_summary ${input_summary}" : ""

    """
    sequenoscope analyze \\
        --input_fastq ${fastqs.join(' ')} \\
        --input_reference ${reference} \\
        --output analyze_results \\
        --sequencing_type ${params.sequencing_type} \\
        ${seq_sum_opt} \\
        ${options}
    """
}

def get_default_analyze(name) {
    def defaults = [
        'minimum_coverage': 1, 'threads': 1, 'minimum_read_length': 15,
        'maximum_read_length': 0, 'trim_front_bp': 0, 'trim_tail_bp': 0,
        'quality_threshold': 15, 'start_time': 0, 'end_time': 100, 'minimap2_kmer': 15
    ]
    return defaults[name]
}

workflow ANALYZE {
    take: ch_to_analyze

    main:
        RUN_ANALYZE(ch_to_analyze)

    emit:
        results = RUN_ANALYZE.out.analysis_results
}