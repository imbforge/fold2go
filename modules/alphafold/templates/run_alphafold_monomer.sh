#!/bin/bash

python /app/alphafold/run_alphafold.py \
    --fasta_paths=${fasta} \
    --output_dir=. \
    --use_precomputed_msas=true \
    --use_gpu_relax=false \
    --max_template_date=2020-05-14 \
    --model_preset=${params.MODEL_PRESET} \
    --data_dir=${params.DATABASE} \
    --bfd_database_path=${params.DATABASE}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --mgnify_database_path=${params.DATABASE}/mgnify/mgy_clusters_2022_05.fa \
    --pdb70_database_path=${params.DATABASE}/pdb70/pdb70 \
    --obsolete_pdbs_path=${params.DATABASE}/pdb_mmcif/obsolete.dat \
    --template_mmcif_dir=${params.DATABASE}/pdb_mmcif/mmcif_files \
    --uniref30_database_path=${params.DATABASE}/uniref30/UniRef30_2021_03 \
    --uniref90_database_path=${params.DATABASE}/uniref90/uniref90.fasta