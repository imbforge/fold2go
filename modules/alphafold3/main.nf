nextflow.preview.types = true

process MSA {
    tag "${meta}"
    label "ssd"

    input:
    (meta, input): Tuple<Map, Set<Path>>

    stage:
    stageAs 'input/*', input

    output:
    json: Path = file("msa/**/*_data.json")

    when:
    params.MSA.enabled

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

    input:
    (meta, json): Tuple<Map, Path>
    cache: Path

    output:
    prediction: Tuple<Map, Path> = tuple(meta, file("predictions/*", type: 'dir'))

    when:
    params.INFERENCE.enabled

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --run_data_pipeline=false \\
        --json_path=${json} \\
        --model_dir=${params.ALPHAFOLD3.MODEL_DIR} \\
        --num_diffusion_samples=${params.ALPHAFOLD3.DIFFUSION_SAMPLES} \\
        --jax_compilation_cache_dir=${cache} \\
        --output_dir=predictions
    """
}