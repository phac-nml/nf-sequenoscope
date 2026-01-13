// modules/filter_ONT.nf
process RUN_FILTER_ONT {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/filter_ONT", 
        mode: 'copy', 
        saveAs: { filename -> file(filename).name }
    
    conda "${moduleDir}/environment.yml"
    // nf-core standard: dynamic selection between native Singularity SIF and Docker image
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequenoscope:1.0.0--pyh7e72e81_1' :
        'quay.io/biocontainers/sequenoscope:1.0.0--pyh7e72e81_1' }"

    input:
    tuple val(meta), path(fastq_inputs), path(input_summary), val(min_ch), val(max_ch)

    output:
    tuple val(meta), path("${meta.id}_filtered/*.fastq"), emit: filtered_results

    script:

    def combined_prefix = params.output_prefix ? "${params.output_prefix}_${meta.id}" : "${meta.id}"

    def final_min = (min_ch != null) ? min_ch : params.minimum_channel // Use provided min_ch in TSV file or default from nextflow.config or CLI
    def final_max = (max_ch != null) ? max_ch : params.maximum_channel // Use provided max_ch in TSV file or default from nextflow.config or CLI

    def channel_options = ""
    if (final_min != null) channel_options += "--minimum_channel ${final_min} "
    if (final_max != null) channel_options += "--maximum_channel ${final_max}"

    log.info "Running Filter_ONT on ${fastq_inputs.join(' ')}"

    """

    sequenoscope filter_ONT \\
        --input_fastq ${fastq_inputs.join(' ')} \\
        --input_summary $input_summary \\
        --output ${meta.id}_filtered \\
        --output_prefix ${combined_prefix} \\
        --minimum_length ${params.minimum_length} \\
        --maximum_length ${params.maximum_length} \\
        --minimum_q_score ${params.minimum_q_score} \\
        --classification ${params.classification} \\
        ${channel_options} \\
        ${params.force ? '--force' : ''}
    """
}


workflow FILTER_ONT {
    take: ch_to_filter
    main:
        RUN_FILTER_ONT(ch_to_filter)
    emit:
        results = RUN_FILTER_ONT.out.filtered_results
}