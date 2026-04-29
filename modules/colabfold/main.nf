nextflow.preview.types = true

process MSA {
    label "ssd"

    input:
    query: List<Path>

    output:
    msa: Path = files("a3m/*.a3m")

    when:
    params.MSA.enabled

    script:
    """
    python ${moduleDir}/resources/usr/bin/combinations.py ${query.join(' ')}

    colabfold_search \\
        --threads ${task.cpus} \\
        combinations.fasta \\
        ${params.COLABFOLD.DATABASE_DIR} \\
        a3m
    """
}

process INFERENCE {
    maxForks 1
    tag "${meta}"
    label "gpu"

    input:
    (meta, msa): Tuple<Map, Path>

    output:
    prediction: Tuple<Map, Path> = tuple(meta, file("predictions", type: 'dir'))
    
    when:
    params.INFERENCE.enabled

    script:
    """
    colabfold_batch \\
        ${msa} \\
        predictions \\
        --num-recycle ${params.COLABFOLD.RECYCLING_STEPS} \\
        --data ${params.COLABFOLD.MODEL_DIR}
    """
}