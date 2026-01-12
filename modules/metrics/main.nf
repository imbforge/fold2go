process METRICS {
    tag "${meta}"
 
    input:
        tuple val(meta), path(prediction)

    output:
        tuple val(meta.model), path("*_metrics.tsv"), emit: metrics

    when:
        params.METRICS.enabled

    script:
        """
        python ${moduleDir}/resources/usr/bin/calculate_metrics.py \\
            --run_name=${workflow.runName} \\
            --predictions=${prediction instanceof List ? '.' : prediction} \\
            --id=${meta.id} \\
            --model_preset=${meta.model}
        """
}
