process ALPHAFOLD {
    label "cuda"
    tag "${meta}"

    when:
        params.ALPHAFOLD.enabled

    input:
        tuple val(meta), path(fasta, stageAs: "chains.fasta"), path(chain_A, stageAs: "chains/msas/A/*"), path(chain_B, stageAs: "chains/msas/B/*")

    output:
        tuple val(meta), path("*.fasta", includeInputs: true), path("chains/*.{pkl,pdb,json}"), emit: prediction

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --fasta_paths=${fasta} \\
            --output_dir=. \\
            --use_precomputed_msas=true \\
            --use_gpu_relax=false \\
            --max_template_date=2020-05-14 \\
            --model_preset=${params.ALPHAFOLD.preset} \\
            --data_dir=${params.ALPHAFOLD.db} \\
            --bfd_database_path=${params.ALPHAFOLD.db}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\
            --mgnify_database_path=${params.ALPHAFOLD.db}/mgnify/mgy_clusters_2022_05.fa \\
            --obsolete_pdbs_path=${params.ALPHAFOLD.db}/pdb_mmcif/obsolete.dat \\
            --pdb_seqres_database_path=${params.ALPHAFOLD.db}/pdb_seqres/pdb_seqres.txt \\
            --template_mmcif_dir=${params.ALPHAFOLD.db}/pdb_mmcif/mmcif_files \\
            --uniprot_database_path=${params.ALPHAFOLD.db}/uniprot/uniprot.fasta \\
            --uniref30_database_path=${params.ALPHAFOLD.db}/uniref30/UniRef30_2021_03 \\
            --uniref90_database_path=${params.ALPHAFOLD.db}/uniref90/uniref90.fasta
        """
}

process MSA {
    tag "${fasta}"

    when:
        params.MSA.enabled

    input:
        path(fasta)
        each(database)

    output:
        tuple val(fasta.baseName), path("*.{a3m,sto}"), emit: msa

    script:
        """
        python ${moduleDir}/resources/usr/bin/run_msa.py \\
            --cores=${task.cpus} \\
            --database=${database} \\
            --fasta_path=${fasta} \\
            --database_root_path=${params.ALPHAFOLD.db}
        """
}
