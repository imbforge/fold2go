process MSA {
    tag "${meta}"

    when:
        params.MSA.enabled

    input:
        tuple val(meta), path(json)

    output:
        tuple val(meta), path("msa/**/*.json"), emit: json

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --norun_inference \\
            --db_dir=${params.DATABASE} \\
            --json_path=${json} \\
            --output_dir=msa
        """
}

process MODEL {
    tag "${meta}"
    label "gpu"

    when:
        params.ALPHAFOLD.enabled

    input:
        tuple val(meta), path(json)

    output:
        tuple val(meta), path("predictions/*.pkl"), emit: prediction

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --norun_data_pipeline \\
            --db_dir=${params.DATABASE} \\
            --model_dir=${params.WEIGHTS}/${params.MODEL_PRESET} \\
            --json_path=${json} \\
            --output_dir=predictions
        """
}