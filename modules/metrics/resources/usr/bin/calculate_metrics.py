#!/usr/bin/env python

import json
import numpy as np

from pathlib import Path
from pandas import DataFrame
from pdb_numpy import Coor
from pdb_numpy.analysis import compute_pdockQ, distance_matrix

def get_contact_pairs(interface: Coor, chains: np.ndarray, cutoff: float = 5.0) -> np.ndarray:
    """
    Get contact pairs at the interface of two chains in a dimeric protein complex.

    Parameters:
    -----------
    interface: Coor
        The interface of the protein complex
    chains: np.ndarray
        The chain identifiers of the two chains
    cutoff: float
        The distance cutoff for contact definition (default: 5.0)
    
    Returns:
    --------
    np.ndarray
        An array of contact pairs at the interface
    """

    receptor = interface.select_atoms(f"chain {chains[0]}")
    contacts = []

    for rec_res in np.unique(receptor.residue):
        ligand = interface.select_atoms(f"chain {chains[1]} and within {cutoff} of residue {rec_res}")
        for lig_res in np.unique(ligand.residue):
            contacts.append([rec_res, lig_res])
    
    return np.array(contacts, dtype=np.int16)

def get_interface_residues(structure: Coor, chains: np.ndarray, cutoff: float = 5.0) -> dict:
    """
    Get the number of residues and atoms at the interface of two chains in a dimeric protein complex.

    Parameters:
    -----------
    structure: Coor
        The protein complex topology
    chains: np.ndarray
        The chain identifiers of the two chains
    cutoff: float
        The distance cutoff for contact definition (default: 5.0)

    Returns:
    --------
    dict
        A dictionary containing the number of residues and atoms at the interface
    """

    # remove hydrogens
    noh = structure.select_atoms("noh")

    # select interface atoms
    interface = noh.select_atoms(
        f"(chain {chains[0]} and within {cutoff} of chain {chains[1]}) or (chain {chains[1]} and within {cutoff} of chain {chains[0]})"
    )

    # select receptor and ligand chains (arbitrary assignment)
    receptor, ligand = interface.select_atoms(f"chain {chains[0]}"), interface.select_atoms(f"chain {chains[1]}")

    # get number of interface residues and atoms
    num_residues, num_atoms = len(get_contact_pairs(interface, chains, cutoff)), np.sum(distance_matrix(receptor.xyz, ligand.xyz) < cutoff)

    return {
        f"num_chain{chains[0]}_intf_res"  : np.unique(receptor.residue).size,
        f"num_chain{chains[1]}_intf_res"  : np.unique(ligand.residue).size,
        'num_res_res_contact'             : num_residues,
        'num_atom_atom_contact'           : num_atoms,
    }

def get_interface_plddt(structure: Coor, chains: np.ndarray, cutoff: float = 5.0) -> dict:
    """
    Get the average pLDDT score at the interface of two chains in a dimeric protein complex.
    Considers all residues with any atom less than 5Å away from any atom of the other chain.
    
    Parameters:
    -----------
    structure: Coor
        The protein complex topology
    chains: np.ndarray
        The chain identifiers of the two chains
    cutoff: float
        The distance cutoff for contact definition (default: 5 Å)

    Returns:
    --------
    dict
        A dictionary containing average pLDDT scores at the interface for each chain and the entire interface
    """

    # remove hydrogens
    noh = structure.select_atoms("noh")

    # select interface atoms
    interface = noh.select_atoms(
        f"(chain {chains[0]} and within {cutoff} of chain {chains[1]}) or (chain {chains[1]} and within {cutoff} of chain {chains[0]})"
    )

    # select receptor and ligand chains (arbitrary assignment)    
    receptor, ligand = interface.select_atoms(f"chain {chains[0]}"), interface.select_atoms(f"chain {chains[1]}")

    # get ipLDDT scores for each chain and the entire interface
    receptor_plddt  = np.fromiter(dict(zip(receptor.residue, receptor.beta)).values(), dtype=np.float64)
    ligand_plddt    = np.fromiter(dict(zip(ligand.residue, ligand.beta)).values(), dtype=np.float64)
    interface_plddt = np.concatenate([receptor_plddt, ligand_plddt])

    return {
        f"chain{chains[0]}_intf_avg_plddt": receptor_plddt.mean(),
        f"chain{chains[1]}_intf_avg_plddt": ligand_plddt.mean(),
        'intf_avg_plddt'                  : interface_plddt.mean()
    }

def get_interface_pae(structure: Coor, chains: np.ndarray, pae: np.ndarray, cutoff: float = 3.5) -> dict:
    """
    Get the median predicted aligned error (iPAE) at the interface of two chains in a dimeric protein complex.
    Following https://github.com/fteufel/alphafold-peptide-receptors/blob/main/qc_metrics.py, the distance threshold to define contact is set to 3.5 Å.

    Parameters:
    -----------
    structure: Coor
        The protein complex topology
    chains: np.ndarray
        The chain identifiers of the two chains
    pae: np.ndarray
        The predicted aligned error for each residue in the protein complex
    cutoff: float
        The distance cutoff for contact definition (default: 3.5 Å)

    Returns:
    --------
    dict
        A dictionary containing the median predicted aligned error at the interface of the protein complex (iPAE)
    """
    
    # select interface atoms
    interface = structure.select_atoms(
        f"(chain {chains[0]} and within {cutoff} of chain {chains[1]}) or (chain {chains[1]} and within {cutoff} of chain {chains[0]})"
    )

    # get contact pairs at the interface
    contacts = get_contact_pairs(interface, chains, cutoff)

    return {
        'iPAE': np.median(pae[np.unique(contacts[:,0]),:][:,np.unique(contacts[:,1])]) if contacts.size else np.nan
    }

def get_pdockq(structure: Coor, chains: np.ndarray, cutoff: float = 8.0) -> dict:
    """
    Get the pDockQ score for a dimeric protein complex.
    Following https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py, the distance threshold to define contact is set to 8 Å.

    Parameters:
    -----------
    structure: Coor
        The protein complex topology
    cutoff: float
        The distance cutoff for contact definition (default: 8.0)
    chain: np.ndarray
        The chain identifiers of the two chains

    Returns:
    --------
    dict
        A dictionary containing the pDockQ score for the protein complex, or mpDockQ if there are more than two chains
    """

    assert chains.size > 1, "At least two chains are required to calculate pDockQ"

    if chains.size == 2:
        return {
            'pDockQ': compute_pdockQ(structure, cutoff=cutoff, L=0.724, x0=152.611, k=0.052, b=0.018).pop()
        }
    else:
        return {
            'mpDockQ': compute_pdockQ(structure, cutoff=cutoff, L=0.728, x0=309.375, k=0.098, b=0.262).pop()
        }

def calculate_af2_metrics(predictions_dir: Path) -> dict:
    """
    Calculate template independent metrics for AlphaFold2 predictions.

    Parameters:
    -----------
    predictions_dir: Path
        The directory containing the predictions

    Returns:
    --------
    dict
        A dictionary containing the calculated metrics for each prediction
    """

    with (predictions_dir / 'ranking_debug.json').open('r') as fin:
        data = json.load(fin)

    ranking = { model:list(data.values())[0].get(model) for model in data['order'] }

    metrics = {}

    for rank, (model, score) in enumerate(ranking.items()):

        coor = Coor(predictions_dir / f"ranked_{rank}.cif")

        chains, lengths = np.unique(coor.select_atoms('protein and name CA').chain, return_counts=True)

        common_metrics = {
            'model_id'               : model,
            'model_rank'             : f"ranked_{rank}",
            'model_confidence'       : float(score),
            **{f"chain{chain}_length": length for chain, length in zip(chains, lengths)},
        }

        if chains.size == 2:
            with (predictions_dir / f"pae_{model}.json").open('r') as fin:
                pae = np.array(json.load(fin)[0]['predicted_aligned_error'], dtype=np.float16)

            metrics[model] = {
                **common_metrics,
                **get_interface_plddt(coor, chains),
                **get_interface_residues(coor, chains),
                **get_interface_pae(coor, chains, pae),
                **get_pdockq(coor, chains)
            }
        else:
            metrics[model] = common_metrics

    return metrics

def calculate_boltz_metrics(predictions_dir: Path) -> dict:
    """
    Calculate template independent metrics for Boltz-1 predictions.

    Parameters:
    -----------
    predictions_dir: Path
        The directory containing the predictions
    
    Returns:
    --------
    dict
        A dictionary containing the calculated metrics for each prediction
    """

    metrics = {}

    for model in predictions_dir.glob('confidence_*.json'):
        with model.open('r') as fin:
            scores = { score: value for score, value in json.load(fin).items() if not isinstance(value, dict) }
        
        model_name = f"{model.parent.name}_model_{model.stem.split('_').pop()}"

        coor = Coor(model.parent / f"{model_name}.cif")

        chains, lengths = np.unique(coor.select_atoms('protein and name CA').chain, return_counts=True)

        common_metrics = {
            'model_id': f"model_{model.stem.split('_').pop()}",
            **{f"chain{chain}_length":length for chain, length in zip(chains, lengths)},
            **scores,
        }

        if chains.size == 2:
            with np.load(model.parent / f"pae_{model_name}.npz") as fin:
                pae = np.array(fin.get('pae'), dtype=np.float16)

            metrics[model_name] = {
                **common_metrics,
                **get_interface_pae(coor, chains, pae),
                **get_interface_residues(coor, chains),
                # **get_interface_plddt(coor, chains), # FIXME: boltz does not populate bfactor
                # **get_pdockq(coor, chains)           # FIXME: boltz does not populate bfactor
            }
        else:
            metrics[model_name] = common_metrics

    return metrics

def calculate_af3_metrics(predictions_dir: Path) -> dict:
    """
    Calculate template independent metrics for AlphaFold3 predictions.

    Parameters:
    -----------
    predictions_dir: Path
        The directory containing the predictions

    Returns:
    --------
    dict
        A dictionary containing the calculated metrics for each prediction
    """

    metrics = {}
    
    for model in predictions_dir.glob('**/seed-*/summary_confidences.json'):

        coor = Coor(model.parent / 'model.cif')

        chains, lengths = np.unique(coor.select_atoms('protein and name CA').chain, return_counts=True)

        with model.open('r') as fin:
            scores = { score: value for score, value in json.load(fin).items() if not isinstance(value, dict | list) }

        common_metrics = {
            'model_id': model.parent.name,
            **{f"chain{chain}_length":length for chain, length in zip(chains, lengths)},
            **scores
        }

        if chains.size == 2:
            with (model.parent / 'confidences.json').open('r') as fin:
                pae = np.array(json.load(fin)['pae'], dtype=np.float16)
            
            metrics[model.parent.name] = {
                    **common_metrics,
                    **get_interface_plddt(coor, chains),
                    **get_interface_residues(coor, chains),
                    **get_interface_pae(coor, chains, pae),
                    **get_pdockq(coor, chains)
            }
        else:
            metrics[model.parent.name] = common_metrics
    
    return metrics

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--predictions', type=Path, dest='predictions')
    parser.add_argument('--id', type=Path, dest='id')
    parser.add_argument('--model_preset', type=str, choices=['alphafold3', 'boltz', 'alphafold2_multimer', 'alphafold2_monomer_ptm', 'alphafold2_monomer', 'alphafold2_monomer_casp14'], dest='model_preset')
    parser.add_argument('--run_name', type=str, dest='run_name')

    args = parser.parse_args()

    meta = {
        'project_name': args.run_name,
        'prediction_name': args.id,
        'model_preset': args.model_preset
    }

    match args.model_preset.split('_')[0]:
        case 'alphafold3':
            metrics = calculate_af3_metrics(args.predictions)
        case 'alphafold2':
            metrics = calculate_af2_metrics(args.predictions)
        case 'boltz':
            metrics = calculate_boltz_metrics(args.predictions)

    DataFrame.from_dict(metrics, orient='index').assign(**meta).round(2).to_csv(f"{args.model_preset}_metrics.tsv", sep='\t', index=False)