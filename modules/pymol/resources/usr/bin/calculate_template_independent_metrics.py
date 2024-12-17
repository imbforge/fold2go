# This script contains generic functions that extract and manipulate metrics and information that can be obtained from predicted models without the need of a template.
# Author: Chop Yan Lee
# Modified by: Patrick Hüther
# pDockQ code source: https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py
# iPAE code source: https://github.com/fteufel/alphafold-peptide-receptors/blob/main/qc_metrics.py

import json, argparse, itertools
import mdtraj as md
import numpy as np
import pandas as pd
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from pymol import cmd
from string import ascii_uppercase

@dataclass
class PredictionFolder:
    """Class that stores prediction folder information"""
    model_preset: str
    prediction_dir: Path
    prediction_name: str
    project_name: str

    @property
    def model_confidences(self) -> dict:
        with open(Path(f'{self.prediction_dir}/ranking_debug.json'), 'r') as fin:
            data = json.load(fin)
        return {model:list(data.values())[0].get(model) for model in data['order']}

    @property
    def chain_lengths(self) -> dict:
        with open(Path(f'{self.prediction_dir}.fasta')) as fin:
            chain_dict, record = defaultdict(int), 0
            for line in fin:
                if line.startswith(">"):
                    chain = f'chain{ascii_uppercase[record]}_length'
                    record += 1
                    continue
                chain_dict[chain] += len(line.strip())
        return chain_dict

    @property
    def model_instances(self) -> list:
        return [
            PredictedModel(model_path=self.prediction_dir,
                           model_id=model,
                           model_rank=f'ranked_{rank}',
                           model_confidence=confidence,
                           model_chains=self.chain_lengths
            ) for rank, (model, confidence) in enumerate(self.model_confidences.items())
        ]

    def write_metrics(self) -> None:
        """Write out the information that has been processed for every predicted model.
        """
        common_metrics = pd.DataFrame({
            'project_name':self.project_name,
            'prediction_name':self.prediction_name,
            'model_preset': f"alphafold2_{self.model_preset}",
            'model_id':self.model_confidences.keys()
        })

        independent_metrics = pd.DataFrame.from_records([
            prediction.model_metrics for prediction in self.model_instances
        ])

        common_metrics.merge(independent_metrics, on='model_id').to_csv('template_indep_info.tsv', sep='\t', index=False)

@dataclass
class PredictedModel:
    """Class that stores predicted model"""
    model_path: Path
    model_id: str
    model_rank: str
    model_confidence: float
    model_chains: dict

    def __post_init__(self):
        self.model_metrics = {
            **self.model_chains,
            'model_id':self.model_id,
            'model_rank':self.model_rank,
            'model_confidence':self.model_confidence
        }

        # calculate additional metrics only for dimers
        if len(self.model_chains) == 2:
            self.chain_coords: dict
            self.chain_plddt: dict

            self.pdb_path: Path = (f'{self.model_path}/{self.model_rank}.pdb')
            self.pae_path: Path = Path(f'{self.model_path}/pae_{self.model_id}.json')

            self.chain_coords, self.chain_plddt = self.parse_pdb()

            self.model_metrics.update({
                **self.calculate_iPAE(),
                **self.calculate_pDockQ(),
                **self.calculate_interface_plddt(),
                **self.calculate_structural_metric()
            })

    def parse_pdb(self) -> tuple:
        """Read a pdb file predicted with AF and rewritten to contain all chains

        Returns:
            self.chain_coords (dict): Dict of chain coordination (x,y,z)
            self.chain_plddt (dict): Dict of chain id as key and plddt array as value
        """
        def _check_chain_id(pdb) -> None:
            """Some models have their chain ids start from B instead of A. As the code requires the chain ids to be consistent (start from chain A), this function checks and rename the chain ids if necessary
            """
            cmd.load(pdb)
            if 'C' in cmd.get_chains():
                # change the chain id into A and B
                cmd.alter('chain B', 'chain="A"')
                cmd.alter('chain C', 'chain="B"')
                # save and overwrite the predicted model
                cmd.save(pdb)
            cmd.reinitialize()

        def _parse_atm_record(line) -> dict:
            """Get the atm record from pdb file

            Returns:
                record (dict): Dict of parsed pdb information from .pdb file
            """
            record = defaultdict()
            record['name'] = line[0:6].strip()
            record['atm_no'] = int(line[6:11])
            record['atm_name'] = line[12:16].strip()
            record['atm_alt'] = line[17]
            record['res_name'] = line[17:20].strip()
            record['chain'] = line[21]
            record['res_no'] = int(line[22:26])
            record['insert'] = line[26].strip()
            record['resid'] = line[22:29]
            record['x'] = float(line[30:38])
            record['y'] = float(line[38:46])
            record['z'] = float(line[46:54])
            record['occ'] = float(line[54:60])
            record['B'] = float(line[60:66])

            return record

        chain_coords, chain_plddt = {}, {}

        _check_chain_id(self.pdb_path)

        with open(self.pdb_path, 'r') as file:
            for line in file:
                if not line.startswith('ATOM'):
                    continue
                record = _parse_atm_record(line)
                #Get CB - CA for GLY
                if record['atm_name']=='CB' or (record['atm_name']=='CA' and record['res_name']=='GLY'):
                    if record['chain'] in [*chain_coords.keys()]:
                        chain_coords[record['chain']].append([record['x'],record['y'],record['z']])
                        chain_plddt[record['chain']].append(record['B'])
                    else:
                        chain_coords[record['chain']] = [[record['x'],record['y'],record['z']]]
                        chain_plddt[record['chain']] = [record['B']]

        #Convert to arrays
        for chain in chain_coords:
            chain_coords[chain] = np.array(chain_coords[chain])
            chain_plddt[chain] = np.array(chain_plddt[chain])

        return (chain_coords, chain_plddt)

    def calculate_pDockQ(self, t=8) -> dict:
        """Wraps the code adapted from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py to calculate pDockQ scores of a given model. Higher score means better.
        pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
        where L= 0.724 x0= 152.611 k= 0.052 and b= 0.018

        Args:
            t (float): distance cutoff (A) to define residues in contact. The authors used 8A between CB atoms or CA for glycine as the cutoff

        Returns:
            pDockQ (float): pDockQ of the model
            pDockQ_PPV (float): Positive predictive value of the given pDockQ
        """
        #Get coords and plddt per chain
        ch1, ch2 = [*self.chain_coords.keys()]
        coords1, coords2 = self.chain_coords[ch1], self.chain_coords[ch2]
        plddt1, plddt2 = self.chain_plddt[ch1], self.chain_plddt[ch2]

        #Calc 2-norm
        mat = np.append(coords1, coords2,axis=0)
        a_min_b = mat[:,np.newaxis,:] - mat[np.newaxis,:,:]
        dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
        l1 = len(coords1)
        contact_dists = dists[:l1,l1:] #upper triangular --> first dim = chain 1
        contacts = np.argwhere(contact_dists<=t)

        if contacts.shape[0] < 1:
            pdockq = 0
            ppv = 0
        else:
            #Get the average interface plDDT
            avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
            #Get the number of interface contacts
            n_if_contacts = contacts.shape[0]
            x = avg_if_plddt*np.log10(n_if_contacts)
            pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611))) + 0.018

            #PPV
            PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192 ,
                0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,
                0.8919553 , 0.88570037, 0.87822061, 0.87116417, 0.86040801,
                0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,
                0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,
                0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,
                0.63555449, 0.55890174])

            pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,
                0.60150931, 0.58313803, 0.5647381 , 0.54122438, 0.52314392,
                0.49659878, 0.4774676 , 0.44661346, 0.42628389, 0.39990988,
                0.38479715, 0.3649393 , 0.34526004, 0.3262589 , 0.31475668,
                0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,
                0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,
                0.06968505, 0.02946438])
            inds = np.argwhere(pdockq_thresholds>=pdockq)
            if len(inds)>0:
                ppv = PPV[inds[-1]][0]
            else:
                ppv = PPV[0]

        return {'pDockQ':pdockq, 'pDockQ_PPV':ppv}

    def calculate_iPAE(self) -> dict:
        """Calculate iPAE using code adapted from https://github.com/fteufel/alphafold-peptide-receptors/blob/main/qc_metrics.py.
        Following the publication, the distance threshold to define contact is set at 0.35nm (3.5A) between CA atoms. Lower score means better.

        Returns:
            iPAE (float): iPAE score of the predicted model
        """
        # Extract plddt and PAE average over binding interface
        model_mdtraj = md.load(self.pdb_path)
        table, _ = model_mdtraj.topology.to_dataframe()
        table = table[(table['name']=='CA')]
        table['residue'] = np.arange(0, len(table))
        
        # receptor (domain) as chainID 0 and ligand (motif) chain as chainID 1 because I always use two chains
        # for prediction and the calling of receptor and ligan is arbitrary
        receptor_res = table[table['chainID'] == 0]['residue']
        ligand_res = table[table['chainID'] == 1]['residue']

        input_to_calc_contacts = [list(product) for product in itertools.product(ligand_res.values,receptor_res.values)]

        contacts, input_to_calc_contacts = md.compute_contacts(model_mdtraj, contacts=input_to_calc_contacts, scheme='closest', periodic=False)
        ligand_res_in_contact = set()
        receptor_res_in_contact = set()

        for i in input_to_calc_contacts[np.where(contacts[0]<0.35)]: # threshold in nm
            ligand_res_in_contact.add(i[0])
            receptor_res_in_contact.add(i[1])
        receptor_res_in_contact = np.fromiter(receptor_res_in_contact, int, len(receptor_res_in_contact))
        ligand_res_in_contact = np.fromiter(ligand_res_in_contact, int, len(ligand_res_in_contact))

        if len(ligand_res_in_contact) > 0:
            with open(self.pae_path, 'rb') as fin:
                pae = np.array(json.load(fin)[0]['predicted_aligned_error'], dtype=np.float16)
            ipae = np.median(pae[receptor_res_in_contact,:][:,ligand_res_in_contact])
        else:
            ipae = 50 # if no residue in contact, impute ipae with large value
        
        return {'iPAE':ipae}

    def calculate_interface_plddt(self) -> dict:
        """Calculate the plddt of every predicted residue using the b-factor of predicted model. Further calculates the average plddt of the residues of each chain at the interface (residues at the interface are defined as the residues with any atom that are less than 5A away from any atom of the other chain). Additionally, create a pymol object that shows the selection of residues at the interface
            
        Returns:
            chainA_intf_plddt (float): average plddt of the residues of chain A that are at the interface
            chainB_intf_plddt (float): average plddt of the residues of chain B that are at the interface
            intf_avg_plddt (float): average plddt of all residues (chain A and B) that are at the interface
            num_chainA_intf_res (int): number of residues of chain A that are at the interface
            num_chainA_intf_res (int): number of residues of chain B that are at the interface
        """
        # load the predicted model
        cmd.load(self.pdb_path)
        # remove hydrogen as they are filled automatically by AlphaFold
        cmd.remove('hydrogens')
        # make a selection of residues in both chain that have at least one atom with less than or equal to 5A from any atom from the other chain
        selection_line = f"({self.model_rank} and chain A within 5A of {self.model_rank} and chain B) or ({self.model_rank} and chain B within 5A of {self.model_rank} and chain A)"
        cmd.select(selection=selection_line, name="residues_less_5A")
        # iterate through the selection by chain to get the b-factor (loaded with plddt by AlphaFold) of the residues
        resi_bfactorA = set()
        resi_bfactorB = set()
        cmd.iterate(f"(residues_less_5A) and chain A","resi_bfactorA.add((resi,b))", space={'resi_bfactorA':resi_bfactorA})
        cmd.iterate(f"(residues_less_5A) and chain B","resi_bfactorB.add((resi,b))", space={'resi_bfactorB':resi_bfactorB})

        # calculate the average plddt of contact residues from each chain
        chainA_intf_plddt = np.array([ele[1] for ele in resi_bfactorA], dtype=np.float16)
        chainB_intf_plddt = np.array([ele[1] for ele in resi_bfactorB], dtype=np.float16)
        intf_plddt = np.concatenate([chainA_intf_plddt, chainB_intf_plddt])

        return {'chainA_intf_avg_plddt':chainA_intf_plddt.mean(), 'chainB_intf_avg_plddt':chainB_intf_plddt.mean(), 'intf_avg_plddt':intf_plddt.mean(), 'num_chainA_intf_res': len(resi_bfactorA), 'num_chainB_intf_res': len(resi_bfactorB)}

    def calculate_structural_metric(self) -> dict:
        """Parse the atoms that are in contact in the predicted model. Contacts are limited to the distance of 5A between two atoms of residues from different chain. Color the predicted model by chain and display the residues in contact as sticks

        Returns:
            num_atom_atom_contact (int): Number of atom-atom contacts
            num_res_res_contact (int): Number of residue-residue contacts
        """
        # make a selection of residues in both chain that have at least one atom with less than or equal to 5A from any atom from the other chain
        selection_line = f"({self.model_rank} and chain A within 5A of {self.model_rank} and chain B) or ({self.model_rank} and chain B within 5A of {self.model_rank} and chain A)"
        cmd.select(selection=selection_line, name="residues_less_5A")
        # iterate through atom indices in one chain that are less than 5A from any atom of the other chain and store them in a list
        chainA_contact_indices = []
        cmd.iterate(f"{self.model_rank} and chain A within 5A of {self.model_rank} and chain B","chainA_contact_indices.append([oneletter,resi,name,index])",space={'chainA_contact_indices':chainA_contact_indices})
        # iterate through the list of indices and find atoms from the other chain that is less than 5A away and calculate distance between them
        unique_resi_pair = set()
        atom_atom_contacts = []
        for oneletterA, resiA, nameA, indexA in chainA_contact_indices:
            chainB_contact_indices = []
            cmd.iterate(f"{self.model_rank} and chain B within 5A from {self.model_rank} and index {indexA}","chainB_contact_indices.append([oneletter,resi,name,index])",space={'chainB_contact_indices':chainB_contact_indices})
            # iterate through the indices from the other chain to calculate their distance with the index from the previous chain
            for oneletterB,resiB,nameB,indexB in chainB_contact_indices:
                dist = cmd.distance("interface_contacts", f"{self.model_rank}`{indexA}",f"{self.model_rank}`{indexB}")
                unique_resi_pair.add((resiA,resiB))
                atom_atom_contact = ['A',oneletterA,resiA,nameA,'B',oneletterB,resiB,nameB,f'{dist:.2f}']
                atom_atom_contacts.append(atom_atom_contact)
        # display the residues in contact as sticks
        cmd.show('cartoon','all')
        cmd.show('sticks','byres (residues_less_5A)')
        cmd.color('green','all and chain A')
        cmd.color('cyan','all and chain B')
        cmd.color('atomic', 'not elem C')
        no_contact_error_message = f'No interface found in {self.model_rank}!'
        try:
            cmd.hide(selection='interface_contact')
        except:
            print(no_contact_error_message)
        cmd.save(f'{self.model_rank}_contacts.pse')
        cmd.reinitialize()

        return {'num_res_res_contact':len(unique_resi_pair), 'num_atom_atom_contact':len(atom_atom_contacts)}

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--model_preset', type=str, choices=['multimer', 'monomer_ptm', 'monomer', 'monomer_casp14'], help='Name of the AlphaFold2 model preset', dest='model_preset')
    parser.add_argument('--prediction_dir', type=Path, help='Path to the prediction directory', dest='prediction_dir')
    parser.add_argument('--prediction_name', type=str, help='Name of the prediction', dest='prediction_name')
    parser.add_argument('--project_name', type=str, help='Name for the project', dest='project_name')

    args = parser.parse_args()

    PredictionFolder(
        model_preset=args.model_preset,
        prediction_dir=args.prediction_dir,
        prediction_name=args.prediction_name,
        project_name=args.project_name
    ).write_metrics()
