#!/bin/bash

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
