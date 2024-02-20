#!/usr/bin/env python

import argparse, shutil

from absl import logging
from alphafold.data.pipeline import run_msa_tool
from alphafold.data.tools import hhblits, jackhmmer

logging.set_verbosity(logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('--database', type=str, choices=['uniref90', 'mgnify', 'uniprot', 'bfd'], help='Name of the database for msa queries')
parser.add_argument('--cores', type=int, help='Number of cpu cores to use')
parser.add_argument('--fasta_path', type=str, help='Path to the input fasta file')
parser.add_argument('--database_root_path', type=str, help='Path to the alphafold database root directory')

args = parser.parse_args()

match args.database:
    case "uniref90":
        database_path=f'{args.database_root_path}/uniref90/uniref90.fasta'
        msa_runner, msa_out_path, max_sto_sequences, msa_format = jackhmmer.Jackhmmer(binary_path=shutil.which('jackhmmer'), database_path=database_path), 'uniref90_hits.sto', 10000, 'sto'
    case "mgnify":
        database_path=f'{args.database_root_path}/mgnify/mgy_clusters_2022_05.fa'
        msa_runner, msa_out_path, max_sto_sequences, msa_format = jackhmmer.Jackhmmer(binary_path=shutil.which('jackhmmer'), database_path=database_path), 'mgnify_hits.sto', 501, 'sto'
    case "uniprot":
        database_path=f'{args.database_root_path}/uniprot/uniprot.fasta'
        msa_runner, msa_out_path, max_sto_sequences, msa_format = jackhmmer.Jackhmmer(binary_path=shutil.which('jackhmmer'), database_path=database_path), 'uniprot_hits.sto', None, 'sto'
    case "bfd":
        databases=[f'{args.database_root_path}/uniref30/UniRef30_2021_03', f'{args.database_root_path}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt']
        msa_runner, msa_out_path, max_sto_sequences, msa_format = hhblits.HHBlits(binary_path=shutil.which('hhblits'), databases=databases), 'bfd_uniref_hits.a3m', None, 'a3m'
    case _:
        parser.error(f"{args.database} is not a valid choice")

msa_runner.n_cpu = args.cores
run_msa_tool(msa_runner=msa_runner, input_fasta_path=args.fasta_path, msa_out_path=msa_out_path, msa_format=msa_format, use_precomputed_msas=False, max_sto_sequences=max_sto_sequences)
