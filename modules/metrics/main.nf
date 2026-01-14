nextflow.preview.types = true

process METRICS {
    tag "${meta}"

    input:
    (meta, prediction): Tuple<Map, Set<Path>>

    output:
    metrics: Tuple<String, Path> = tuple(meta.model, file("*_metrics.tsv"))

    when:
    params.METRICS.enabled

    script:
    """
    python ${moduleDir}/resources/usr/bin/calculate_metrics.py \\
        --run_name=${workflow.runName} \\
        --predictions=${prediction.size() > 1 ? '.' : prediction.pop()} \\
        --id=${meta.id} \\
        --model_preset=${meta.model}
    """
}
