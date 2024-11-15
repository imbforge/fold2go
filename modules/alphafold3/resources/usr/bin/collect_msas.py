#!/usr/bin/env python

import argparse
import json
from hashlib import md5
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--json_path', type=Path, required=True)
parser.add_argument('--msa_dir', type=Path, default=Path.cwd() / 'msa')
args = parser.parse_args()

def _get_msa(entity, realm, msa_dir=args.msa_dir):
    seqhash = md5(entity[realm]['sequence'].encode('utf-8')).hexdigest()
    try:
        with open(f'{msa_dir}/{seqhash}_data.json') as fin:
            record = json.load(fin)['sequences'][0][realm]
    except FileNotFoundError:
        return { 'unpairedMsa': "" }
    match (realm):
        case 'rna':
            return { 'unpairedMsa': record['unpairedMsa'] }
        case 'protein':
            return { key: record[key] for key in ['unpairedMsa', 'templates'] }

with open(args.json_path) as fin:
    jobdef = json.load(fin)

for i, entity in enumerate(jobdef['sequences']):
    if (realm := list(entity).pop()) in ['protein','rna']:
        jobdef['sequences'][i][realm] |= _get_msa(entity, realm)
    else:
        continue

with open(f"{jobdef['name']}_data.json", 'w') as fout:
    json.dump(jobdef, fout, indent=2)
