process MSA {
    tag "${meta}"

    when:
        params.MSA.enabled

    input:
        tuple val(meta), path(input)

    output:
        path("msa/**/*_data.json"), emit: json

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --run_inference=false \\
            --${input.extension == "json" ? "json_path" : "input_dir"}=${input} \\
            --db_dir=${params.DATABASE_DIR}/${params.ALPHAFOLD_VERSION} \\
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
            --model_dir=${params.MODEL_DIR} \\
            --jax_compilation_cache_dir=${workDir} \\
            --output_dir=predictions
        """

}