process METRICS {
    tag "${meta}"
 
    when:
        params.METRICS.enabled

    input:
        tuple val(meta), path(prediction)

    output:
        tuple val(meta.model), path("*_metrics.tsv"), emit: metrics

    script:
        """
        python ${moduleDir}/resources/usr/bin/calculate_metrics.py \\
            --run_name=${workflow.runName} \\
            --predictions=${prediction instanceof List ? '.' : prediction} \\
            --id=${meta.id} \\
            --model_preset=${meta.model}
        """
}
