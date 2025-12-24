// modules/filter_ONT.nf
process RUN_FILTER_ONT {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/filter_ONT", 
        mode: 'copy', 
        saveAs: { filename -> file(filename).name }
    conda "bioconda::sequenoscope=1.0.0"

    input:
    tuple val(meta), path(fastq_inputs), path(input_summary), val(reference), val(min_ch), val(max_ch)

    output:
    tuple val(meta), path("filter_results/*.fastq"), val(reference), emit: filtered_results

    script:
    def flags = [
        'minimum_channel', 'maximum_channel', 'minimum_duration', 
        'maximum_duration', 'minimum_start_time', 'maximum_start_time',
        'minimum_q_score', 'maximum_q_score', 'minimum_length', 'maximum_length',
        'output_prefix', 'classification'
    ]

    def options = flags.collect { flag ->
        if (params[flag] != null && params[flag] != get_default(flag)) {
            return "--${flag} ${params[flag]}"
        }
    }.findAll().join(' ')

    if (params.summarize) options += " --summarize"
    if (params.force)     options += " --force"

    """
    sequenoscope filter_ONT \\
        --input_fastq ${fastq_inputs.join(' ')} \\
        --input_summary $input_summary \\
        --output filter_results \\
        $options
    """
}

def get_default(name) {
    def defaults = [
        'classification': 'all', 'minimum_channel': 1, 'maximum_channel': 512, 
        'minimum_duration': 0, 'maximum_duration': 100, 'minimum_start_time': 0, 
        'maximum_start_time': 259200, 'minimum_q_score': 0, 'maximum_q_score': 100, 
        'minimum_length': 0, 'maximum_length': 50000, 'output_prefix': 'sample'
    ]
    return defaults[name]
}

workflow FILTER_ONT {
    take: ch_to_filter
    main:
        RUN_FILTER_ONT(ch_to_filter)
    emit:
        results = RUN_FILTER_ONT.out.filtered_results
}