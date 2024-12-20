process INFERENCE_MONOMER {
    tag "${meta}"
    label "gpu"

    when:
        params.INFERENCE.enabled

    input:
        tuple val(meta), path(fasta), path(chain, stageAs: "chain/msas/*")

    output:
        tuple val(meta), path("chain/*.{pdb,pkl,cif,json}"), emit: prediction

    script:
        template 'run_alphafold_monomer.sh'
}

process INFERENCE_MULTIMER {
    tag "${meta}"
    label "gpu"

    when:
        params.INFERENCE.enabled

    input:
        tuple val(meta),
              path(fasta,  stageAs: "chains.fasta"   ),
              path(chainA, stageAs: "chains/msas/A/*"),
              path(chainB, stageAs: "chains/msas/B/*"),
              path(chainC, stageAs: "chains/msas/C/*"),
              path(chainD, stageAs: "chains/msas/D/*"),
              path(chainE, stageAs: "chains/msas/E/*"),
              path(chainF, stageAs: "chains/msas/F/*"),
              path(chainG, stageAs: "chains/msas/G/*"),
              path(chainH, stageAs: "chains/msas/H/*")

    output:
        tuple val(meta), path("chains/*.{pdb,pkl,cif,json}"), emit: prediction

    script:
        template 'run_alphafold_multimer.sh'
}

process MSA {
    tag "${record.id}:${database}"
    label "ssd"

    when:
        params.MSA.enabled

    input:
        tuple val(meta), val(record)
        each(database)

    output:
        tuple val(record.id), path("msas/*/*.{a3m,sto}"), emit: msa

    script:
        def chain = meta.find { it.value == record.id }.key
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
