nextflow.enable.types = true

process INFERENCE {
    tag "${fasta.baseName}"
    label "gpu"

    input:
    record(
        chains: Map,
        A     : Set<Path>,
        B     : Set<Path>,
        C     : Set<Path?>,
        D     : Set<Path?>,
        E     : Set<Path?>,
        F     : Set<Path?>,
        G     : Set<Path?>,
        H     : Set<Path?>,
        fasta : Path,
        model : String
    )

    stage:
    stageAs A, "chains/msas/A/*"
    stageAs B, "chains/msas/B/*"
    stageAs C, "chains/msas/C/*"
    stageAs D, "chains/msas/D/*"
    stageAs E, "chains/msas/E/*"
    stageAs F, "chains/msas/F/*"
    stageAs G, "chains/msas/G/*"
    stageAs H, "chains/msas/H/*"

    output:
    record(
        id        : fasta.baseName,
        prediction: file("chains", type: 'dir'),
        model     : model
    )

    script:
    """
    mv ${fasta} chains/chains.fasta

    python /app/alphafold/run_alphafold.py \\
        --fasta_paths=chains/chains.fasta \\
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

    rm -rf chains/msas
    """
}

process MSA {
    tag "${id}:${db}"
    label "ssd"

    input:
    record(
        chains: Map,
        id    : String,
        fasta : Path,
        db    : String,
        model : String
    )

    output:
    record(
        chains: chains,
        id    : id,
        msa   : file("msas/${id}/*.{a3m,sto}"),
        model : model
    )

    script:
    """
    python ${moduleDir}/resources/usr/bin/run_msa.py \\
        --cores=${task.cpus} \\
        --database=${db} \\
        --fasta_path=${fasta} \\
        --data_dir=${params.ALPHAFOLD2.DATABASE_DIR} \\
        --out_path=msas/${id}
    """
}
