nextflow.enable.types = true

process MSA {
    label "ssd"

    input:
    record(
        model: String,
        query: Set<Path>
    )

    stage:
    stageAs query, 'input/*'

    output:
    record(
        model: model,
        msa  : files("**/*_data.json")
    )

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --run_inference=false \\
        --input_dir=input \\
        --db_dir=${params.ALPHAFOLD3.DATABASE_DIR} \\
        --output_dir=msa
    """
}

process INFERENCE {
    tag "${id}"
    label "gpu"

    input:
    record(
        id   : String,
        msa  : Path,
        model: String
    )

    output:
    record(
        id        : id,
        prediction: file("predictions/${id}", type: 'dir'),
        model     : model
    )

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --run_data_pipeline=false \\
        --json_path=${msa} \\
        --model_dir=${params.ALPHAFOLD3.MODEL_DIR} \\
        --num_diffusion_samples=${params.ALPHAFOLD3.DIFFUSION_SAMPLES} \\
        --jax_compilation_cache_dir=${workflow.workDir} \\
        --output_dir=predictions
    """
}