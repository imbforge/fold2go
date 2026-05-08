nextflow.enable.types = true

process MSA {
    label "ssd"

    input:
    record(
        model: String,
        query: Set<Path>
    )

    output:
    record(
        model: model,
        msa  : files("msa/*.a3m"),
        json : files("msa/*.json", optional: true)
    )

    script:
    """
    python ${moduleDir}/resources/usr/bin/combinations.py ${query.join(' ')}

    colabfold_search \\
        --threads ${task.cpus} \\
        combinations.fasta \\
        ${params.COLABFOLD.DATABASE_DIR} \\
        ${model == 'alphafold3' ? '--af3-json' : ''} \\
        msa
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
        prediction: file("predictions", type: 'dir'),
        model     : model
    )

    script:
    """
    colabfold_batch \\
        ${msa} \\
        predictions \\
        --num-recycle ${params.COLABFOLD.RECYCLING_STEPS} \\
        --data ${params.COLABFOLD.MODEL_DIR}
    """
}