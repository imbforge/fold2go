nextflow.preview.types = true

process INFERENCE_MONOMER {
    tag "${meta}"
    label "gpu"

    input:
    (meta, fasta, chain): Tuple<Map, Path, List<Path>>

    stage:
    stageAs "chain/msas/*", chain

    output:
    prediction: Tuple<Map, Set<Path>> = tuple(meta, files("chain/*.{pdb,pkl,cif,json}"))

    when:
    params.INFERENCE.enabled

    script:
    """
    grep -A1 ${meta['A']} ${fasta} > chain.fasta

    python /app/alphafold/run_alphafold.py \\
        --fasta_paths=chain.fasta \\
        --output_dir=. \\
        --use_precomputed_msas=true \\
        --use_gpu_relax=false \\
        --max_template_date=2020-05-14 \\
        --model_preset=${params.ALPHAFOLD2.MODEL_PRESET} \\
        --data_dir=${params.ALPHAFOLD2.DATABASE_DIR} \\
        --bfd_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\
        --mgnify_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/mgnify/mgy_clusters_2022_05.fa \\
        --pdb70_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/pdb70/pdb70 \\
        --obsolete_pdbs_path=${params.ALPHAFOLD2.DATABASE_DIR}/pdb_mmcif/obsolete.dat \\
        --template_mmcif_dir=${params.ALPHAFOLD2.DATABASE_DIR}/pdb_mmcif/mmcif_files \\
        --uniref30_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/uniref30/UniRef30_2023_02 \\
        --uniref90_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/uniref90/uniref90.fasta
    """
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
    prediction: Tuple<Map, Set<Path>> = tuple(meta, files("chains/*.{pdb,pkl,cif,json}"))

    when:
    params.INFERENCE.enabled

    script:
    """
    python /app/alphafold/run_alphafold.py \\
        --fasta_paths=${fasta} \\
        --output_dir=. \\
        --use_precomputed_msas=true \\
        --use_gpu_relax=false \\
        --max_template_date=2020-05-14 \\
        --model_preset=multimer \\
        --num_multimer_predictions_per_model=${params.ALPHAFOLD2.PREDICTIONS_PER_MODEL} \\
        --data_dir=${params.ALPHAFOLD2.DATABASE_DIR} \\
        --bfd_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\
        --mgnify_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/mgnify/mgy_clusters_2022_05.fa \\
        --pdb_seqres_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/pdb_seqres/pdb_seqres.txt \\
        --obsolete_pdbs_path=${params.ALPHAFOLD2.DATABASE_DIR}/pdb_mmcif/obsolete.dat \\
        --template_mmcif_dir=${params.ALPHAFOLD2.DATABASE_DIR}/pdb_mmcif/mmcif_files \\
        --uniprot_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/uniprot/uniprot.fasta \\
        --uniref30_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/uniref30/UniRef30_2023_02 \\
        --uniref90_database_path=${params.ALPHAFOLD2.DATABASE_DIR}/uniref90/uniref90.fasta
    """
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
