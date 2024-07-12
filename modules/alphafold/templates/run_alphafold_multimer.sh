#!/bin/bash

python /app/alphafold/run_alphafold.py \\
    --fasta_paths=${fasta} \\
    --output_dir=. \\
    --use_precomputed_msas=true \\
    --use_gpu_relax=false \\
    --max_template_date=2020-05-14 \\
    --model_preset=multimer \\
    --num_multimer_predictions_per_model=${params.PREDICTIONS_PER_MODEL} \\
    --data_dir=${params.DATABASE} \\
    --bfd_database_path=${params.DATABASE}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\
    --mgnify_database_path=${params.DATABASE}/mgnify/mgy_clusters_2022_05.fa \\
    --pdb_seqres_database_path=${params.DATABASE}/pdb_seqres/pdb_seqres.txt \\
    --obsolete_pdbs_path=${params.DATABASE}/pdb_mmcif/obsolete.dat \\
    --template_mmcif_dir=${params.DATABASE}/pdb_mmcif/mmcif_files \\
    --uniprot_database_path=${params.DATABASE}/uniprot/uniprot.fasta \\
    --uniref30_database_path=${params.DATABASE}/uniref30/UniRef30_2023_02 \\
    --uniref90_database_path=${params.DATABASE}/uniref90/uniref90.fasta
