process RUN_PLOT {
    tag "${meta.id}"
    // Following your 'analyze' pattern: results go into a specific plot folder
    publishDir "${params.output}/${meta.id}/plots", 
        mode: 'copy', 
        saveAs: { filename -> file(filename).name }
    
    conda "bioconda::sequenoscope=1.0"

    input:
    // We pass the meta object so we know which sample we are plotting
    tuple val(meta), path('test_input_dir'), path('ctrl_input_dir')

    output:
    tuple val(meta), path("plot_results/*"), emit: plot_results

    script:
    // These names now match your nextflow.config keys exactly
    def plot_flags = [
        'output_prefix', 
        'violin_data_percent', 
        'time_bin_unit'
    ]

    def options = plot_flags.collect { flag ->
        if (params[flag] != null && params[flag] != get_default_plot(flag)) {
            return "--${flag} ${params[flag]}"
        }
    }.findAll().join(' ')

    // Handle boolean flags (flags that don't take a value)
    if (params.adaptive_sampling) options += " --adaptive_sampling"
    if (params.force)             options += " --force"

    """
    sequenoscope plot \\
        --test_dir test_input_dir \\
        --control_dir ctrl_input_dir \\
        --output_dir plot_results \\
        ${options}
    """
}

def get_default_plot(name) {
    def defaults = [
        'output_prefix': 'sample', 
        'violin_data_percent': 0.1, 
        'time_bin_unit': 'minutes'
    ]
    return defaults[name]
}


workflow PLOT {
    take: ch_plot_input
    main:
        RUN_PLOT(ch_plot_input)
    emit:
        results = RUN_PLOT.out.plot_results
}