process MSA {
    tag "${meta}"
    label "ssd"

    when:
        params.MSA.enabled

    input:
        tuple val(meta), path(input, stageAs: 'input/*')

    output:
        path("msa/**/*_data.json"), emit: json

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --run_inference=false \\
            --${input instanceof List ? "input_dir=input" : "json_path=" << input} \\
            --db_dir=${params.ALPHAFOLD3.DATABASE_DIR} \\
            --output_dir=msa
        """
}

process INFERENCE {
    tag "${meta}"
    label "gpu"

    when:
        params.INFERENCE.enabled

    input:
        tuple val(meta), path(json)

    output:
        tuple val(meta), path("predictions/*", type: 'dir'), emit: prediction

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --run_data_pipeline=false \\
            --json_path=${json} \\
            --model_dir=${params.ALPHAFOLD3.MODEL_DIR} \\
            --jax_compilation_cache_dir=${workDir} \\
            --output_dir=predictions
        """

}