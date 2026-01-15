nextflow.preview.types = true

process MSA {
    tag "${meta}"
    label "ssd"

    input:
    (meta, input): Tuple<Map, Set<Path>>

    stage:
    stageAs 'input/*', input

    output:
    msa: Tuple<Map, Path> = tuple([ id: meta.id, model: 'alphafold3' ], file("**/**/*_data.json"))

    when:
    params.MSA.enabled

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --run_inference=false \\
        --${input.size() > 1 ? "input_dir=input" : "json_path=" << input.pop()} \\
        --db_dir=${params.ALPHAFOLD3.DATABASE_DIR} \\
        --output_dir=msa
    """
}

process INFERENCE {
    tag "${meta}"
    label "gpu"

    input:
    (meta, json): Tuple<Map, Path>

    output:
    prediction: Tuple<Map, Path> = tuple(meta, file("predictions/${meta.id}", type: 'dir'))

    when:
    params.INFERENCE.enabled

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --run_data_pipeline=false \\
        --json_path=${json} \\
        --model_dir=${params.ALPHAFOLD3.MODEL_DIR} \\
        --num_diffusion_samples=${params.ALPHAFOLD3.DIFFUSION_SAMPLES} \\
        --jax_compilation_cache_dir=${workflow.workDir} \\
        --output_dir=predictions
    """
}