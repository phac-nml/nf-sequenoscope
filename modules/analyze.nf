process RUN_ANALYZE {
    tag "${meta.id}"

    publishDir "${params.output}/${meta.id}/", mode: 'copy', overwrite: true

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/sequenoscope:1.0.0--pyh7e72e81_1" // Nextflow automatically prefixes this with 'docker://' when using Singularity

    input:
    tuple val(meta), path(fastqs), path(input_summary), path(reference), val(min_ch), val(max_ch)

    output:
    tuple val(meta), path("${meta.id}_analyze_${meta.group}results"), emit: analysis_results

    script:
    def combined_prefix = params.output_prefix ? "${params.output_prefix}_${meta.id}" : "${meta.id}"
    def seq_type = (fastqs instanceof List && fastqs.size() == 2) ? 'PE' : params.sequencing_type

    // FIX: Correctly split fastq files for paired-end mode
    def fastq1 = ""
    def fastq2 = ""

    if (fastqs instanceof List && fastqs.size() == 2) {
        fastq1 = fastqs[0]
        fastq2 = fastqs[1]
    } else {
        fastq1 = fastqs instanceof List ? fastqs[0] : fastqs
    }
    // Handle the sequencing summary and type explicitly as they are often required
    def seq_sum_opt = (input_summary && !input_summary.empty()) ? "--sequencing_summary ${input_summary}" : ""
    log.info "Running ANALYZE on ID: ${meta.id} | Mode: ${seq_type} | Files: ${fastqs} | Ref: ${reference}"

    """
    sequenoscope analyze \\
        --input_fastq ${fastq1} ${fastq2} \\
        --input_reference ${reference} \\
        --output ${meta.id}_analyze_${meta.group}results \\
        --sequencing_type ${seq_type} \\
        --output_prefix ${combined_prefix} \\
        --sequencing_type ${seq_type} \\
        --threads ${params.threads} \\
        --minimum_coverage ${params.minimum_coverage} \\
        --minimum_read_length ${params.minimum_read_length} \\
        --maximum_read_length ${params.maximum_read_length} \\
        --trim_front_bp ${params.trim_front_bp} \\
        --trim_tail_bp ${params.trim_tail_bp} \\
        --quality_threshold ${params.quality_threshold} \\
        ${seq_sum_opt} \\
        ${params.force ? '--force' : '--force'}
    """
}



workflow ANALYZE {
    take: ch_to_analyze

    main:
        RUN_ANALYZE(ch_to_analyze)

    emit:
        results = RUN_ANALYZE.out.analysis_results
}