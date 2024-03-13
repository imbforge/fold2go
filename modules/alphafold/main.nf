process ALPHAFOLD {
    label "cuda"
    tag "${meta}"

    when:
        params.ALPHAFOLD.enabled

    input:
        tuple val(meta), path(chain_A, stageAs: "chains/msas/A/*"), path(chain_B, stageAs: "chains/msas/B/*"), path(fasta, stageAs: "chains.fasta")

    output:
        tuple val(meta), path("chains/*.{pkl,pdb,json}"), emit: prediction

    script:
        """
        python /app/alphafold/run_alphafold.py \\
            --fasta_paths=${fasta} \\
            --output_dir=. \\
            --use_precomputed_msas=true \\
            --use_gpu_relax=false \\
            --max_template_date=2020-05-14 \\
            --model_preset=${params.MODEL_PRESET} \\
            --num_multimer_predictions_per_model=${params.PREDICTIONS_PER_MODEL} \\
            --data_dir=${params.DATABASE} \\
            --bfd_database_path=${params.DATABASE}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\
            --mgnify_database_path=${params.DATABASE}/mgnify/mgy_clusters_2022_05.fa \\
            --obsolete_pdbs_path=${params.DATABASE}/pdb_mmcif/obsolete.dat \\
            --pdb_seqres_database_path=${params.DATABASE}/pdb_seqres/pdb_seqres.txt \\
            --template_mmcif_dir=${params.DATABASE}/pdb_mmcif/mmcif_files \\
            --uniprot_database_path=${params.DATABASE}/uniprot/uniprot.fasta \\
            --uniref30_database_path=${params.DATABASE}/uniref30/UniRef30_2021_03 \\
            --uniref90_database_path=${params.DATABASE}/uniref90/uniref90.fasta
        """
}

process MSA {
    tag "${record.id}"

    when:
        params.MSA.enabled

    input:
        tuple val(meta), val(record)
        each(database)

    output:
        tuple val(record.id), path("*.{a3m,sto}"), emit: msa

    script:
        """
        cat << EOF > ${record.id}.fasta
        >chain_${meta.find { it.value == record.id }.key}
        ${record.seqString}
        EOF

        python ${moduleDir}/resources/usr/bin/run_msa.py \\
            --cores=${task.cpus} \\
            --database=${database} \\
            --fasta_path=${record.id}.fasta \\
            --database_root_path=${params.DATABASE}
        """
}
