# This script contains generic functions that extract and manipulate metrics and information that can be obtained from predicted models without the need of a template.
# Author: Chop Yan Lee
# pDockQ code source: https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py
# iPAE code source: https://github.com/fteufel/alphafold-peptide-receptors/blob/main/qc_metrics.py

from pymol import cmd
import numpy as np
import pandas as pd
import json, os, pickle, argparse, sys, itertools
from collections import defaultdict
import mdtraj as md

class Prediction_folder:
    """Class that stores prediction folder information"""
    def __init__(self,prediction_folder,num_model=5,project_name=None):
        """Initialize an instance of Prediction

        Args:
            prediction_folder (str): absolute path to the prediction folder
        """
        self.prediction_folder = prediction_folder
        self.num_model = num_model
        self.path_to_prediction_folder = os.path.split(self.prediction_folder)[0]
        self.prediction_name = os.path.split(self.prediction_folder)[1]
        self.rank_to_model = {}
        self.model_confidences = {}
        self.fasta_sequence_dict = {'A':'','B':''}
        # instantiate the amount of Predicted_model according to the number of models given as argument, otherwise 5
        self.model_instances = {}
        if project_name is not None:
            self.project_name = project_name
        # need an attribute to annotate if a prediction folder has been successfully predicted without internal AlphaFold error
        self.predicted = True

    def parse_ranking_debug_file(self):
        """Read the ranking_debug_file and save relevant information into attribute of self
        """
        if not os.path.exists(os.path.join(self.prediction_folder,'ranking_debug.json')):
            self.predicted = False
            return
        else:
            with open(os.path.join(self.prediction_folder,'ranking_debug.json'), 'r') as f:
                data = json.load(f)
            self.rank_to_model = {f'ranked_{i}':model for i, model in enumerate(data.get("order"))}
            sorted_model_confidence = sorted(data.get("iptm+ptm").values(),reverse=True)
            self.model_confidences = {f'ranked_{i}':float(confidence) for i, confidence in enumerate(sorted_model_confidence)}
        
    def parse_prediction_fasta_file(self):
        """Read the fasta file of the prediction to retrieve information on chain and sequence identity
        """
        fasta_path = f'{self.prediction_folder}.fasta'
        with open(fasta_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip() != '']
        chain_id = 0
        for line in lines:
            if line[0] == '>':
                chain = list(self.fasta_sequence_dict)[chain_id]
                chain_id += 1
                continue
            self.fasta_sequence_dict[chain] += line

    def instantiate_predicted_model(self):
        """Initialize the amount of Predicted_model instance according to the number of model specified and save it in the dict self.model_instances
        """
        self.model_instances = {f'ranked_{i}':Predicted_model(f'ranked_{i}') for i in range(self.num_model)}

    def assign_model_info(self):
        """Assign information stored in the prediction folder to their corresponding predicted model
        """
        for model_id, model_inst in self.model_instances.items():
            model_inst.model_confidence = self.model_confidences.get(model_id)
            model_inst.multimer_model = self.rank_to_model.get(model_id)
            model_inst.path_to_model = self.prediction_folder

    def process_all_models(self):
        """Use the instances of Predicted_model and run the wrapper function Predicted_model.get_model_independent_metrics function on themselves
        """
        self.parse_ranking_debug_file()
        self.parse_prediction_fasta_file()
        if self.predicted:
            self.instantiate_predicted_model()
            self.assign_model_info()
            for model_id, model_inst in self.model_instances.items():
                model_inst.get_model_independent_metrics()

    def write_out_calculated_metrics(self,project_name=None):
        """Write out the information that has been processed for every predicted model. Check if template_indep_info.tsv already exists in the same level as all prediction folder, if it does, read in the file as pd.Dataframe and append new info into it, otherwise create a new one

        Args:
            project_name (str): a project name given to model contacts dataframe as a key identifier

        Returns:
            template_indep_info.tsv: A tsv file with the calculated template independent metrics
        """
        metrics_out_path = os.path.join(self.path_to_prediction_folder,'template_indep_info.tsv')
        # prepare the metrics dataframe
        if self.project_name is not None:
            metrics_columns_dtype = {'project_name':str}
        else:
            metrics_columns_dtype = {}
        # include the rest of the columns
        metrics_columns_dtype.update({
            'prediction_name':str,
            'chainA_length':int,
            'chainB_length':int,
            'model_id':str,
            'model_confidence':float,
            'chainA_intf_avg_plddt':float,
            'chainB_intf_avg_plddt':float,
            'intf_avg_plddt':float,
            'pDockQ':float,
            'iPAE':float,
            'num_chainA_intf_res':int,
            'num_chainB_intf_res':int,
            'num_res_res_contact':int,
            'num_atom_atom_contact':int
            })
        # check if template_indep_info.tsv already exists
        if os.path.exists(metrics_out_path):
            metrics_df = pd.read_csv(metrics_out_path,sep='\t',index_col=0)
            metrics_df.reset_index(drop=True,inplace=True)
        else:
            metrics_df = pd.DataFrame(columns=metrics_columns_dtype.keys())
            metrics_df = metrics_df.astype(dtype=metrics_columns_dtype)

        if self.project_name is not None:
            common_info = [self.project_name]
        else:
            common_info = []
        common_info = common_info + [self.prediction_name,len(self.fasta_sequence_dict.get('A')), len(self.fasta_sequence_dict.get('B'))]
        # check if the prediction folder has been predicted successfully without internal error from AlphaFold
        if not self.predicted:
            row = common_info + ['Prediction failed'] + [None]*10
            metrics_df.loc[len(metrics_df)] = row
        else:
            # insert metric info in a row-wise manner
            for model_id, model_inst in self.model_instances.items():
                row = common_info + [model_id]
                for column in list(metrics_columns_dtype)[5:]:
                    row.append(model_inst.__dict__.get(column))
                metrics_df.loc[len(metrics_df)] = row
        metrics_df.to_csv(metrics_out_path,sep='\t')
        print(f'Calculated metrics saved in {metrics_out_path}!')

    def write_out_contacts(self):
        """Write out atom-atom contacts into a tsv file
        
        Returns:
            {self.run_id}_all_model_contacts.tsv: A tsv file with the calculated information of atom-atom contact in all predicted models
        """
        if self.predicted:
            contact_out_path = os.path.join(self.path_to_prediction_folder,'all_model_contacts.tsv')
            # prepare the model contact dataframe
            contact_columns_dtype = {
                'project_name':str,
                'prediction_name':str,
                'model_id':str,
                'chain1':str,
                'res1':str,
                'resID1':int,
                'atom1':str,
                'chain2':str,
                'res2':str,
                'resID2':int,
                'atom2':str,
                'atom_distance':float
                }
            # check if all_model_contacts.tsv already exists
            if os.path.exists(contact_out_path):
                contact_df = pd.read_csv(contact_out_path,sep='\t',index_col=0)
                contact_df.reset_index(drop=True,inplace=True)
            else:
                contact_df = pd.DataFrame(columns=contact_columns_dtype.keys())
                contact_df = contact_df.astype(dtype=contact_columns_dtype)
            
            # insert the contact info in a row-wise manner
            common_info = [self.project_name,self.prediction_name]
            for model_id, model_inst in self.model_instances.items():
                for atom_atom_contact in model_inst.atom_atom_contacts:
                    row = common_info + [model_id] + atom_atom_contact
                    contact_df.loc[len(contact_df)] = row
            contact_df.to_csv(contact_out_path,sep='\t')
            print(f'Calculated atom-atom contacts saved in {contact_out_path}!')


class Predicted_model:
    """Class that stores predicted model"""
    def __init__(self,predicted_model):
        """Initialize an instance of Predicted_model
        
        Args:
            predicted_model (str): name of the predicted model like ranked_0
        """
        self.predicted_model = predicted_model
        self.path_to_model = None
        self.multimer_model = None
        self.chain_coords = None
        self.chain_plddt = None
        self.pickle_data = None
        self.model_confidence = None
        self.chainA_intf_avg_plddt = None
        self.chainB_intf_avg_plddt = None
        self.intf_avg_plddt = np.nan
        self.pDockQ = np.nan
        self.PPV = np.nan
        self.iPAE = np.nan
        self.num_chainA_intf_res = 0
        self.num_chainB_intf_res = 0
        self.num_res_res_contact = 0
        self.atom_atom_contacts = []
        self.num_atom_atom_contact = 0

    def read_pickle(self):
        """Read in the pickle file of multimer model

        Returns:
            self.pickle_data (dict): Pickle data of multimer model
        """
        multimer_model_pickle = os.path.join(self.path_to_model,f'result_{self.multimer_model}.pkl')
        with open(multimer_model_pickle, 'rb') as f:
            self.pickle_data = pickle.load(f)

    def check_chain_id(self):
        """Some models have their chain ids start from B instead of A. As the code requires the chain ids to be consistent (start from chain A), this function checks and rename the chain ids if necessary
        """
        model_path = os.path.join(self.path_to_model,f'{self.predicted_model}.pdb')
        cmd.load(model_path)
        chains = cmd.get_chains(f'{self.predicted_model}')
        if 'C' in chains:
            # change the chain id into a temporary arbitrary name
            cmd.alter(f'{self.predicted_model} and chain B', 'chain="tempA"')
            cmd.sort()
            cmd.alter(f'{self.predicted_model} and chain C', 'chain="tempB"')
            cmd.sort()
            # change the chain id into A and B
            cmd.alter(f'{self.predicted_model} and chain tempA', 'chain="A"')
            cmd.sort()
            cmd.alter(f'{self.predicted_model} and chain tempB', 'chain="B"')
            cmd.sort()
            # save and overwrite the predicted model
            cmd.save(model_path,f'{self.predicted_model}')
        cmd.reinitialize()

    def parse_atm_record(self,line):
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
    
    def read_pdb(self):
        """Read a pdb file predicted with AF and rewritten to conatin all chains

        Returns:
            self.chain_coords (dict): Dict of chain coordination (x,y,z)
            self.chain_plddt (dict): Dict of chain id as key and plddt array as value
        """

        chain_coords, chain_plddt = {}, {}
        model_path = os.path.join(self.path_to_model,f'{self.predicted_model}.pdb')

        with open(model_path, 'r') as file:
            for line in file:
                if not line.startswith('ATOM'):
                    continue
                record = self.parse_atm_record(line)
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

        self.chain_coords = chain_coords
        self.chain_plddt = chain_plddt

    def calc_pdockq(self, t=8):
        """Calculate the pDockQ scores 
        pdockQ = L / (1 + np.exp(-k*(x-x0)))+b
        where L= 0.724 x0= 152.611 k= 0.052 and b= 0.018

        Args:
            t (float): distance cutoff (A) to define residues in contact. The authors used 8A between CB atoms or CA for glycine as the cutoff

        Returns:
            self.pDockQ (float): pDockQ of the model
            self.PPV (float): Positive predictive value of the given pDockQ
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

        if contacts.shape[0]<1:
            pdockq=0
            ppv=0
        else:
            #Get the average interface plDDT
            avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
            #Get the number of interface contacts
            n_if_contacts = contacts.shape[0]
            x = avg_if_plddt*np.log10(n_if_contacts)
            pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

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

        self.pDockQ = pdockq
        self.PPV = ppv

    def calculate_pDockQ(self):
        """Wraps the code adapted from https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py to calculate pDockQ of a given model. Higher score means better.
        """
        self.read_pdb()
        self.calc_pdockq()

    def calculate_iPAE(self):
        """Calculate iPAE using code adapted from https://github.com/fteufel/alphafold-peptide-receptors/blob/main/qc_metrics.py. Following the publication, the distance threshold to define contact is set at 0.35nm (3.5A) between CA atoms. Lower score means better.

        Returns:
            self.iPAE (float): iPAE score of the predicted model
        """
        multimer_model_pickle = os.path.join(self.path_to_model,f'result_{self.multimer_model}.pkl')
        model_path = os.path.join(self.path_to_model,f'{self.predicted_model}.pdb')

        df = pd.DataFrame({})

        prediction = pd.read_pickle(multimer_model_pickle)

        # Extract plddt and PAE average over binding interface
        model_mdtraj = md.load(model_path)
        table, _ = model_mdtraj.topology.to_dataframe()
        table = table[(table['name']=='CA')]
        table['residue'] = np.arange(0, len(table))
        
        # receptor (domain) as chainID 0 and ligand (motif) chain as chainID 1 because I always use two chains
        # for prediction and the calling of receptor and ligan is arbitrary
        receptor_res = table[table['chainID'] == 0]['residue']
        ligand_res = table[table['chainID'] == 1]['residue']

        input_to_calc_contacts = [list(product) for product in itertools.product(ligand_res.values,receptor_res.values)]

        contacts, input_to_calc_contacts = md.compute_contacts(model_mdtraj, contacts=input_to_calc_contacts,scheme='closest', periodic=False)
        ligand_res_in_contact = []
        receptor_res_in_contact = []

        for i in input_to_calc_contacts[np.where(contacts[0]<0.35)]: # threshold in nm
            ligand_res_in_contact.append(i[0])
            receptor_res_in_contact.append(i[1])
        receptor_res_in_contact, receptor_res_counts = np.unique(np.asarray(receptor_res_in_contact),return_counts=True)
        ligand_res_in_contact, ligand_res_counts = np.unique(np.asarray(ligand_res_in_contact), return_counts=True)

        if len(ligand_res_in_contact) > 0:
            ipae = np.median(prediction['predicted_aligned_error'][receptor_res_in_contact,:][:,ligand_res_in_contact])
        else:
            ipae = 50 # if no residue in contact, impute ipae with large value
        
        self.iPAE = ipae

    def parse_ptm_iptm(self):
        """Parse the ptm and iptm of a predicted model by using the pickle file of the multimer model where the ptm and iptm can be found
        
        Returns:
            ptm (float): the parsed ptm, saved as attribute of self
            iptm (float): the parsed iptm, saved as attribute of self
        """
        self.ptm = float(self.pickle_data['ptm'])
        self.iptm = float(self.pickle_data['iptm'])

    def calculate_interface_plddt(self):
        """Calculate the plddt of every predicted residue using the b-factor of predicted model. Further calculates the average plddt of the residues of each chain at the interface (residues at the interface are defined as the residues with any atom that are less than 5A away from any atom of the other chain). Additionally, create a pymol object that shows the selection of residues at the interface
            
        Returns:
            chainA_intf_plddt (float): average plddt of the residues of chain A that are at the interface, saved as attribute of self
            chainB_intf_plddt (float): average plddt of the residues of chain B that are at the interface, saved as attribute of self
            intf_avg_plddt (float): average plddt of all residues (chain A and B) that are at the interface, saved as attribute of self
        """
        model_path = os.path.join(self.path_to_model,f'{self.predicted_model}.pdb')
        # load the predicted model
        cmd.load(model_path)
        # remove hydrogen as they are filled automatically by AlphaFold
        cmd.extract('h_atoms', 'hydrogens', source_state=1, target_state=1)
        cmd.delete('h_atoms')
        cmd.remove('hydrogen')
        cmd.sort()
        # make a selection of residues in both chain that have at least one atom with less than or equal to 5A from any atom from the other chain
        selection_line = f"({self.predicted_model} and chain A within 5A of {self.predicted_model} and chain B) or ({self.predicted_model} and chain B within 5A of {self.predicted_model} and chain A)"
        cmd.select(selection=selection_line, name="residues_less_5A")
        # iterate through the selection by chain to get the b-factor (loaded with plddt by AlphaFold) of the residues
        resi_bfactorA = set()
        resi_bfactorB = set()
        cmd.iterate(f"(residues_less_5A) and chain A","resi_bfactorA.add((resi,b))",space={'resi_bfactorA':resi_bfactorA})
        cmd.iterate(f"(residues_less_5A) and chain B","resi_bfactorB.add((resi,b))",space={'resi_bfactorB':resi_bfactorB})
        # add number of interface residue to the attribute
        self.num_chainA_intf_res = len(resi_bfactorA)
        self.num_chainB_intf_res = len(resi_bfactorB)
        # calculate the average plddt of contact residues from each chain
        if any(resi_bfactorA):
            self.chainA_intf_avg_plddt = np.mean([float(ele[1]) for ele in resi_bfactorA])
            self.chainB_intf_avg_plddt = np.mean([float(ele[1]) for ele in resi_bfactorB])
            self.intf_avg_plddt = np.mean([float(ele[1]) for ele in resi_bfactorA] + [float(ele[1]) for ele in resi_bfactorB])

    def calculate_structural_metric(self):
        """Parse the atoms that are in contact in the predicted model. Contacts are limited to the distance of 5A between two atoms of residues from different chain. Color the predicted model by chain and display the residues in contact as sticks
            
        Returns:
            self.atom_atom_contacts (list of list): A list of atom-atom contacts stored as list [chain_idA, residueA, residue_indexA, atom_nameA, chain_idB, residueB, residue_indexB, atom_nameB, distance]
            self.num_res_res_contact (int): Number of residue-residue contacts
        """
        # make a selection of residues in both chain that have at least one atom with less than or equal to 5A from any atom from the other chain
        selection_line = f"({self.predicted_model} and chain A within 5A of {self.predicted_model} and chain B) or ({self.predicted_model} and chain B within 5A of {self.predicted_model} and chain A)"
        cmd.select(selection=selection_line, name="residues_less_5A")
        # iterate through atom indices in one chain that are less than 5A from any atom of the other chain and store them in a list
        chainA_contact_indices = []
        cmd.iterate(f"{self.predicted_model} and chain A within 5A of {self.predicted_model} and chain B","chainA_contact_indices.append([oneletter,resi,name,index])",space={'chainA_contact_indices':chainA_contact_indices})
        # iterate through the list of indices and find atoms from the other chain that is less than 5A away and calculate distance between them
        unique_resi_pair = set()
        atom_atom_contacts = []
        for oneletterA,resiA,nameA,indexA in chainA_contact_indices:
            chainB_contact_indices = []
            cmd.iterate(f"{self.predicted_model} and chain B within 5A from {self.predicted_model} and index {indexA}","chainB_contact_indices.append([oneletter,resi,name,index])",space={'chainB_contact_indices':chainB_contact_indices})
            # iterate through the indices from the other chain to calculate their distance with the index from the previous chain
            for oneletterB,resiB,nameB,indexB in chainB_contact_indices:
                dist = cmd.distance("interface_contacts",f"{self.predicted_model}`{indexA}",f"{self.predicted_model}`{indexB}")
                unique_resi_pair.add((resiA,resiB))
                atom_atom_contact = ['A',oneletterA,resiA,nameA,'B',oneletterB,resiB,nameB,f'{dist:.2f}']
                atom_atom_contacts.append(atom_atom_contact)
        # display the residues in contact as sticks
        cmd.show('cartoon','all')
        cmd.show('sticks','byres (residues_less_5A)')
        cmd.color('green','all and chain A')
        cmd.color('cyan','all and chain B')
        cmd.color('atomic', 'not elem C')
        no_contact_error_message = f'No interface found in {self.predicted_model}!'
        try:
            cmd.hide(selection='interface_contact')
        except:
            print(no_contact_error_message)
        self.num_res_res_contact = len(unique_resi_pair)
        self.atom_atom_contacts = atom_atom_contacts
        self.num_atom_atom_contact = len(atom_atom_contacts)
        cmd.save(f'{os.path.join(self.path_to_model,self.predicted_model)}_contacts.pse')
        print(f'PyMol session of calculated contacts saved in {os.path.join(self.path_to_model,self.predicted_model)}_contacts.pse!')
        cmd.reinitialize()

    def get_model_independent_metrics(self):
        """Wraps all the functions together to process a predicted model

        Returns:
            None
        """
        # self.parse_ptm_iptm() # skipped for now because the pickle file has JAX dependency and I am not sure what to do with it
        try:
            self.read_pickle()
            self.calculate_iPAE()
        except:
            print("Error loading the pickle file, skipping the calculation of iPAE.")
        # if 'multimer_v2' in self.multimer_model:
        #     if os.path.exists(os.path.join(self.path_to_model,f'result_{self.multimer_model}.pkl')):
        #         self.read_pickle()
        #         self.calculate_iPAE()
        self.check_chain_id()
        self.calculate_pDockQ()
        self.calculate_interface_plddt()
        self.calculate_structural_metric()
        print(f'{os.path.join(self.path_to_model,self.predicted_model)} processed!')

def main():
    """Parse arguments and wraps all functions into main for executing the program in such a way that it can handle multiple run ids given to it
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-run_ids', type=str, help='Run IDs for metrics calculation', dest='run_ids')
    parser.add_argument('-path_to_run', type=str, help='Either provide a path to a folder where multiple runs of AlphaFold predictions are contained and specify the run_ids to be processed or use -path_to_prediction to specify a folder that you want to process, include "/" at the end', dest='path_to_run')
    parser.add_argument('-path_to_prediction', type=str, help='Path to the prediction folder "/" at the end', dest='path_to_prediction')
    parser.add_argument('-project_name', type=str, help='Optional name for the project', dest='project_name')
    parser.add_argument('-skip_write_out_contacts', action='store_true', help='Exclude writing out atom-atom contacts found in predicted models', dest='skip_write_out_contacts')
    args = parser.parse_args()
    run_ids = vars(args)['run_ids']
    path_to_run = vars(args)['path_to_run']
    path_to_prediction = vars(args)['path_to_prediction']
    project_name = vars(args)['project_name']
    skip_contacts = vars(args)['skip_write_out_contacts']

    # a list to contains already processed files
    calculated_files = []

    # check which argument, -path_to_run or -path_to_prediction, is provided
    if (path_to_run is None) and (path_to_prediction is None):
        print('Please provide either -path_to_run or -path_to_prediction and try again!')
        sys.exit()
    elif path_to_prediction is not None:
        if os.path.exists(f'{path_to_prediction}template_indep_info.tsv'):
            temp = pd.read_csv(f'{path_to_prediction}template_indep_info.tsv',sep='\t',index_col=0)
            calculated_files = temp['prediction_name'].unique()
        for file in os.listdir(path_to_prediction):
            if file in calculated_files:
                # print(file)
                continue
            file_abs = os.path.join(path_to_prediction,file)
            if os.path.isdir(file_abs):
                folder = Prediction_folder(file_abs,num_model=5,project_name=project_name)
                folder.process_all_models()
                folder.write_out_calculated_metrics()
                if skip_contacts:
                    continue
                folder.write_out_contacts()
    else:
        for run_id in run_ids.split(','):
            if os.path.exists(f'{path_to_run}run{run_id}/template_indep_info.tsv'):
                temp = pd.read_csv(f'{path_to_run}run{run_id}/template_indep_info.tsv',sep='\t',index_col=0)
                calculated_files = temp['prediction_name'].unique()
            run_path = f'{path_to_run}run{run_id}'
            for file in os.listdir(run_path):
                if file in calculated_files:
                    # print(file)
                    continue
                file_abs = os.path.join(run_path,file)
                if os.path.isdir(file_abs):
                    folder = Prediction_folder(file_abs,num_model=5,project_name=project_name)
                    folder.process_all_models()
                    folder.write_out_calculated_metrics()
                    if skip_contacts:
                        continue
                    folder.write_out_contacts()

if __name__ == '__main__':
    main()

# python3 ~/AlphaFold_benchmark/scripts/calculate_template_independent_metrics.py -run_ids 41 -path_to_run /Volumes/imb-luckgr/projects/dmi_predictor/DMI_AF2_PRS/ -project_name AlphaFold_benchmark
# python3 ~/AlphaFold_benchmark/scripts/calculate_template_independent_metrics.py -run_ids 41 -path_to_run /Volumes/imb-luckgr/projects/dmi_predictor/DMI_AF2_PRS/ -project_name AlphaFold_benchmark -skip_write_out_contacts