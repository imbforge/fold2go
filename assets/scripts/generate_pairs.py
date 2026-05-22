#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

def parse_fasta(file_path: Path):
	def __parse_uniprot_id(header: str):
		# try to extract uniprot_id from fasta headers
		# return original header otherwise
		if "|" in header:
			parts = header.split("|")
			if (acc := parts[1].strip()):
				return acc
		return header

	records = {}
	with file_path.open("r", encoding="utf-8") as fasta:
		for line in fasta:
			if line.startswith('>'):
				header = __parse_uniprot_id(line.strip()[1:])
				records[header] = ''
			elif header is not None:
				records[header] += line.rstrip("\n")
			else:
				raise ValueError(f"{file_path} is not a valid FASTA file!")
	return records

def build_payload(id: str, seq1: str, seq2: str):
	return {
		"name": id,
		"sequences": [
			{"protein": {"id": "A", "sequence": seq1}},
			{"protein": {"id": "B", "sequence": seq2}},
		],
		"modelSeeds": [1],
		"dialect": "alphafold3",
		"version": 1,
	}

def main():
	parser = argparse.ArgumentParser(
		description=(
			"Read two FASTA files and build cartesian product of their records. Supports fasta format for colabfold (mmseqs2) and json format for af3 (jackhmmer)"
		)
	)
	parser.add_argument("fasta1", type=Path, help="First FASTA file")
	parser.add_argument("fasta2", type=Path, help="Second FASTA file")
	parser.add_argument(
		"format",
		type=str,
		help="Choose fasta or json output",
		default='fasta',
		const='fasta',
		nargs='?',
		choices=('fasta', 'json')
	)
	parser.add_argument(
		"-o",
		"--output-dir",
		type=Path,
		default=Path("."),
		help="Directory where JSON files are written (default: current directory)",
	)
	args = parser.parse_args()

	records1 = parse_fasta(args.fasta1)
	records2 = parse_fasta(args.fasta2)

	args.output_dir.mkdir(parents=True, exist_ok=True)

	# get all combinations and store them in dictionary
	combinations = {
		f"{acc1}_vs_{acc2}": {"seq1": seq1, "seq2": seq2}
		for acc1, seq1 in records1.items()
		for acc2, seq2 in records2.items()
	}

	match args.format:
		case 'json':
			for id, seqs in combinations.items():
				output_file = args.output_dir / f"{id}.json"
				payload = build_payload(id=id, seq1=seqs['seq1'], seq2=seqs['seq2'])
				with output_file.open("w", encoding="utf-8") as handle:
					json.dump(payload, handle, indent=4)
					handle.write("\n")
		case 'fasta':
			output_file = args.output_dir / f"{list(records1.keys())[0]}_vs_many.fasta"
			with output_file.open("w", encoding="utf-8") as handle: 
				for id, seqs in combinations.items():
					handle.write(f">{id}\n{seqs['seq1']}:{seqs['seq2']}\n")

if __name__ == "__main__":
	main()

