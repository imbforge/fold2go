nextflow.enable.types = true

process METRICS {
    tag "${id}"

    input:
    record(
        id: String,
        prediction: Path,
        model: String
    )

    output:
    record(id: id, metrics: file("*_metrics.tsv"), model: model)

    script:
    """
    python ${moduleDir}/resources/usr/bin/calculate_metrics.py \\
        --run_name=${workflow.runName} \\
        --predictions=${prediction} \\
        --id=${id} \\
        --model_preset=${model}
    """
}
