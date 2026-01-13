nextflow.preview.types = true

process INFERENCE_MONOMER {
    tag "${meta}"
    label "gpu"

    input:
    (meta, _fasta, chain): Tuple<Map, Path, List<Path>>

    stage:
    stageAs "chain/msas/*", chain

    output:
    prediction: Tuple<Map, Path> = tuple(meta, file("chain", type: 'dir'))

    when:
    params.INFERENCE.enabled

    script:
    template('run_alphafold_monomer.sh')
}

process INFERENCE_MULTIMER {
    tag "${meta}"
    label "gpu"

    input:
    (meta, fasta, chainA, chainB, chainC, chainD, chainE, chainF, chainG, chainH): Tuple<Map, Path, List<Path>, List<Path>, List<Path?>, List<Path?>, List<Path?>, List<Path?>, List<Path?>, List<Path?>>

    stage:
    stageAs "chains.fasta", fasta
    stageAs "chains/msas/A/*", chainA
    stageAs "chains/msas/B/*", chainB
    stageAs "chains/msas/C/*", chainC
    stageAs "chains/msas/D/*", chainD
    stageAs "chains/msas/E/*", chainE
    stageAs "chains/msas/F/*", chainF
    stageAs "chains/msas/G/*", chainG
    stageAs "chains/msas/H/*", chainH

    output:
    prediction: Tuple<Map, Path>  = tuple(meta, file("chains", type: 'dir'))

    when:
    params.INFERENCE.enabled

    script:
    template('run_alphafold_multimer.sh')
}

process MSA {
    tag "${record.id}:${database}"
    label "ssd"

    input:
    (meta, record, database): Tuple<Map, Map, String>

    output:
    msa: Tuple<Map, Path> = tuple(record.id, file("msas/*/*.{a3m,sto}"))

    when:
    params.MSA.enabled

    script:
    def chain = meta.find { m -> m.value == record.id }.key
    """
        cat << EOF > '${record.id}.fasta'
        >chain_${chain}
        ${record.seqString}
        EOF

        python ${moduleDir}/resources/usr/bin/run_msa.py \\
            --cores=${task.cpus} \\
            --database=${database} \\
            --fasta_path=${record.id}.fasta \\
            --data_dir=${params.ALPHAFOLD2.DATABASE_DIR} \\
            --out_path=msas/${record.id}
        """
}
