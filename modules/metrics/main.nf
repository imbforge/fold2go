nextflow.preview.types = true

process METRICS {
    tag "${meta}"

    input:
    (meta, prediction): Tuple<Map, Path>

    output:
    metrics: Tuple<String, Path> = tuple(meta, file("*_metrics.tsv"))

    when:
    params.METRICS.enabled

    script:
    """
    python ${moduleDir}/resources/usr/bin/calculate_metrics.py \\
        --run_name=${workflow.runName} \\
        --predictions=${prediction} \\
        --id=${meta.id} \\
        --model_preset=${meta.model}
    """
}
