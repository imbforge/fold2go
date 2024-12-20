process METRICS {
    tag "${meta}"
 
    when:
        params.METRICS.enabled

    input:
        tuple val(meta), path(prediction)

    output:
        tuple val(meta.model), path("*_metrics.tsv"), emit: metrics

    script:
        println prediction
        """
        python ${moduleDir}/resources/usr/bin/calculate_metrics.py \\
            --run_name=${workflow.runName} \\
            --predictions=${prediction} \\
            --id=${meta.id} \\
            --model_preset=${meta.model}
        """
}
