process RUN_PLOT {
    tag "${meta.id}"

    publishDir "${params.output}/${meta.id}", 
        mode: 'copy', 
        saveAs: { filename -> file(filename).name }
    
    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/sequenoscope:1.0.0--pyh7e72e81_1" // Nextflow automatically prefixes this with 'docker://' when using Singularity

    input:
    tuple val(meta), path(test_dir), path(ctrl_dir)

    output:
    tuple val(meta), path("plot_results/*"), emit: plot_results

    script:

    // Handle boolean flags (flags that don't take a value)
    def combined_prefix = params.output_prefix ? "${params.output_prefix}_${meta.id}" : "${meta.id}"
    log.info "Running PLOT on test_dir: '${test_dir}' and ctrl_dir: '${ctrl_dir}'"

    """
    sequenoscope plot \\
        --test_dir ${test_dir} \\
        --control_dir ${ctrl_dir} \\
        --output_dir plot_results \\
        --output_prefix ${combined_prefix} \\
        --violin_data_percent ${params.violin_data_percent} \\
        --time_bin_unit ${params.time_bin_unit} \\
        ${params.adaptive_sampling ? '--adaptive_sampling' : ''} \\
        ${params.force ? '--force' : '--force'}
    """
}


workflow PLOT {
    take: ch_plot_input
    main:
        RUN_PLOT(ch_plot_input)
    emit:
        results = RUN_PLOT.out.plot_results
}