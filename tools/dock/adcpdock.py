import os
import shutil
import tempfile
import subprocess
import dataclasses as dc
from typing import List, Optional
from Bio import PDB
from Bio.PDB import Model as PDBModel
from ppflow.datasets.constants import *
from tools.relax.rosetta_packing import side_chain_packing
from Bio.PDB import PDBParser, PDBIO
from tqdm import tqdm
import numpy as np
import pandas as pd


parser = PDBParser(QUIET=True)


import abc
from typing import List


FilePath = str


class DockingEngine(abc.ABC):

    @abc.abstractmethod
    def __enter__(self):
        pass

    @abc.abstractmethod
    def __exit__(self, typ, value, traceback):
        pass

    @abc.abstractmethod
    def set_receptor(self, pdb_path: FilePath):
        pass

    @abc.abstractmethod
    def set_ligand(self, pdb_path: FilePath):
        pass

    @abc.abstractmethod
    def dock(self) -> List[FilePath]:
        pass

def fix_docked_pdb(pdb_path):
    fixed = []
    with open(pdb_path, 'r') as f:
        for ln in f.readlines():
            if (ln.startswith('ATOM') or ln.startswith('HETATM')) and len(ln) == 56:
                fixed.append( ln[:-1] + ' 1.00  0.00              \n' )
            else:
                fixed.append(ln)
    with open(pdb_path, 'w') as f:
        f.write(''.join(fixed))


def extract_energy_info(save_dir, output_file='energy_summary.txt'):
    energy_file = os.path.join(save_dir, "energy_record.txt")
    if not os.path.exists(energy_file):
        print("Energy record file not found.")
        return

    metadata = []
    energy_data = []
    table_started = False
    headers = []

    with open(energy_file, 'r') as f:
        for line in f:
            if line.startswith("mode |"):
                table_started = True
                headers = line.strip().split()
                continue
            if table_started:
                if line.strip() and not line.startswith("-----+"):
                    energy_data.append(line.split())
            elif line.strip():
                metadata.append(line.strip())

    # Create a DataFrame from the energy data
    df = pd.DataFrame(energy_data)

    # Assign headers to the DataFrame
    if len(headers) < len(df.columns):
        df.columns = headers + [f'Unnamed_{i}' for i in range(len(df.columns) - len(headers))]
    else:
        df.columns = headers[:len(df.columns)]

    # Convert numeric columns to float
    numeric_columns = ['mode', 'affinity', 'clust.', 'ref.', 'rmsd', 'energy', 'best']
    for col in df.columns:
        if any(nc in col for nc in numeric_columns):
            try:
                df[col] = df[col].apply(lambda x: pd.to_numeric(x, errors='coerce'))
            except Exception as e:
                print(f"Error converting column {col} to numeric: {e}")
                print(f"Column contents: {df[col].tolist()[:5]}...")  # Print first 5 elements for debugging

    # Sort by affinity (assuming lower is better)
    df = df.sort_values('affinity')

    # Calculate summary statistics for numerical columns
    summary = df.select_dtypes(include=[np.number]).describe()

    # Write metadata, table, and summary to the output file
    with open(os.path.join(save_dir, output_file), 'w') as f:
        # Write metadata
        f.write("Metadata:\n")
        for line in metadata:
            f.write(f"{line}\n")
        f.write("\n")

        # Write energy table
        f.write("Energy Table:\n")
        f.write(df.to_string(index=False))
        f.write("\n\n")

        # Write summary statistics
        f.write("Summary Statistics:\n")
        f.write(summary.to_string())

    print(f"Energy information saved to {os.path.join(save_dir, output_file)}")





class ADCPDock(DockingEngine):

    def __init__(
            self,
            save_dir,
            reduce_bin = None,
            adcp_bin = None,
            prepare_ligand_bin = None,
            prepare_receptor_bin = None,
            agfr_bin = None
        ):
        super().__init__()
        
        # Get the current working directory
        current_dir = os.getcwd()
        
        # Set default paths relative to the current directory if not provided
        if reduce_bin is None:
            reduce_bin = os.path.join(current_dir, 'bin', 'ADFRsuite_x86_64Linux_1.0', 'bin', 'reduce')
        reduce_het_dict = os.path.join(current_dir, 'bin/ADFRsuite_x86_64Linux_1.0/bin/reduce_wwPDB_het_dict.txt')
        if adcp_bin is None:
            adcp_bin = os.path.join(current_dir, 'bin', 'ADFRsuite_x86_64Linux_1.0', 'bin', 'adcp')
        if prepare_ligand_bin is None:
            prepare_ligand_bin = os.path.join(current_dir, 'bin', 'ADFRsuite_x86_64Linux_1.0', 'bin', 'prepare_ligand')
        if prepare_receptor_bin is None:
            prepare_receptor_bin = os.path.join(current_dir, 'bin', 'ADFRsuite_x86_64Linux_1.0', 'bin', 'prepare_receptor')
        if agfr_bin is None:
            agfr_bin = os.path.join(current_dir, 'bin', 'ADFRsuite_x86_64Linux_1.0', 'bin', 'agfr')

        # Use os.path.realpath to resolve any symbolic links
        self.adcp_bin = os.path.realpath(adcp_bin)
        self.reduce_bin = os.path.realpath(reduce_bin)
        self.reduce_het_dict = os.path.realpath(reduce_het_dict)
        self.prepare_ligand_bin = os.path.realpath(prepare_ligand_bin)
        self.prepare_receptor_bin = os.path.realpath(prepare_receptor_bin)
        self.agfr_bin = os.path.realpath(agfr_bin)

        self.tmpdir = tempfile.TemporaryDirectory()
        self.save_dir = save_dir

        self._has_receptor = False
        self._has_ligand = False
        self._ligand_packed = False
        self._receptor_packed = False

    def __enter__(self):

        return self
    
    def __exit__(self, typ, value, traceback):

        self.tmpdir.cleanup()

    def set_receptor(self, pdb_path):

        shutil.copyfile(pdb_path, os.path.join(self.tmpdir.name, 'receptor.pdb'))
        self._receptor_path = os.path.join(self.tmpdir.name,'receptor.pdb')

        self._has_receptor = True

    def set_ligand(self, pdb_path):

        shutil.copyfile(pdb_path, os.path.join(self.tmpdir.name, 'ligand.pdb'))
        self._ligand_path = os.path.join(self.tmpdir.name, 'ligand.pdb')
        self._has_ligand = True

    def side_chain_packing(self, type):
        out_path = side_chain_packing(pdb_file=os.path.join(self.tmpdir.name, '{}.pdb'.format(type)),
                                      output_file=os.path.join(self.tmpdir.name, '{}_packed.pdb'.format(type)))
        
        if type == 'ligand':
            self._ligand_path = out_path
            self._ligand_packed = True
        elif type == 'receptor':
            self._receptor_path = out_path
            self._receptor_packed = True

        else:
            raise ValueError('No such type')

        return out_path


    def _dump_complex_pdb(self, save_name, n_save):
        
        shutil.copyfile(os.path.join(self.tmpdir.name,  "rl_redocking_{}.pdb".format(n_save)),
                        os.path.join(self.save_dir, "{}_{}.pdb".format(save_name, n_save)))
        
        shutil.copyfile(os.path.join(self.tmpdir.name, "terminal_record.txt"),
                        os.path.join(self.save_dir, f"{save_name}_energy_record.txt"))
        
        energies = []
        
        with open(os.path.join(self.tmpdir.name, "terminal_record.txt"), 'r') as f:

            data = f.read().splitlines()
            
            energs = []

            for i, line in enumerate(data):

                if line[0] == "-":

                    energs = data[i+1 : -1]

            for e in energs:

                energies.append(float(e.split()[1]))

        if len(energies) == 0:
            print('Failed!')
            return np.array([np.nan])

        return np.array(energies)

        

    def load_pep_seq_and_struc(self):
        structure_ligand = parser.get_structure(self._ligand_path, self._ligand_path)[0]

        for chain in structure_ligand:
            seq = []
            struc = []
            for residue in chain:
                seq.append(resindex_to_ressymb[AA(residue.get_resname())])
                atom_names = ['N', 'CA', 'C']
                
                for idx, atom_name in enumerate(atom_names):
                    if atom_name == '': continue
                    if atom_name in residue:
                        struc.append(np.array(residue[atom_name].get_coord().tolist()))   

        seq = ''.join(seq)
        seq = seq.lower()

        return seq, np.array(struc)

    def _reduce_receptor(self, receptor_h_path):
        reduce_receptor_list = [self.reduce_bin, os.path.basename(self._receptor_path), '>',
                        os.path.basename(receptor_h_path)]
        cmdline_cd = 'cd {}'.format(self.tmpdir.name)
        cmdline = cmdline_cd + ' && ' + ' '.join(reduce_receptor_list)
        os.system(cmdline)

    def _reduce_ligand(self, ligand_h_path):
        reduce_ligand_list = [self.reduce_bin, os.path.basename(self._ligand_path), '>',
                        os.path.basename(ligand_h_path)]
        cmdline_cd = 'cd {}'.format(self.tmpdir.name)
        cmdline = cmdline_cd + ' && ' + ' '.join(reduce_ligand_list)
        os.system(cmdline)

    def _bio_reduce_ligand(self, ligand_h_path):
        pdb_parser = PDB.PDBParser(QUIET=True)
        structure = pdb_parser.get_structure('protein', self._ligand_path)

        hydrogen_pdb = PDB.PDBIO()
        hydrogen_pdb.set_structure(structure)

        hydrogen_pdb.save(ligand_h_path, select=None)

    def auto_dock_box1(self, lig_struc):
        docking_center = lig_struc.mean(0)
        box_size = np.abs(lig_struc - docking_center[None,:]).max(0) * 2.0  # Increased multiplier
        padding = 10.0  # Add some padding
        box_size += padding
        return docking_center, box_size
    
    
    def auto_dock_box(self, receptor_structure, target_chain, target_residues):
        # Parse the receptor structure
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("receptor", receptor_structure)

        
        # Extract coordinates of target residues
        target_coords = []
        found_residues = set()
        chain = structure[0][target_chain]  # Get first model, specified chain

        for residue in chain:
                if residue.id[1] in target_residues:
                    found_residues.add(residue.id[1])
                    for atom in residue:
                        if atom.name == 'CA':  # Use alpha carbon as reference
                            target_coords.append(atom.coord)
            
        if not target_coords:
            print(f"Warning: No target residues found in chain {target_chain}! Falling back to ligand-based box")
            lig_seq, lig_struc = self.load_pep_seq_and_struc()
            return self.auto_dock_box1(lig_struc)
        
        print(f"Found residues in chain {target_chain}: {found_residues}")
        
        target_coords = np.array(target_coords)
        docking_center = np.mean(target_coords, axis=0)
        box_size = np.max(np.abs(target_coords - docking_center), axis=0) * 1.5
        padding = 20.0
        box_size += padding
        
        print(f"Docking box center: {docking_center}")
        print(f"Docking box size: {box_size}")
        
        return docking_center, box_size

    def dock(self, save_name, n_save, n_search=20, n_steps=1000000, auto_box=False):

        if not (self._has_receptor and self._has_ligand):

            raise ValueError('Missing receptor or peptide.')

        lig_seq, lig_struc = self.load_pep_seq_and_struc()
        
        if auto_box:
            # Specify which chain and residues to use
            
            target_chain = 'C'  # Change this to your desired chain
            target_residues = [59, 98, 112]  # Your target residues
            
            try:
                docking_center, box_size = self.auto_dock_box(self._receptor_path, target_chain, target_residues)
            except Exception as e:
                print(f"Error in auto_dock_box: {e}")
                print("Falling back to ligand-based box")
                docking_center, box_size = self.auto_dock_box1(lig_struc)


        receptor_h_path = self._receptor_path.split('/')[-1].split('.')[0] + 'H.pdb'
        receptor_h_path = os.path.join(self.tmpdir.name, receptor_h_path)

        receptor_pt_path = self._receptor_path.split('/')[-1].split('.')[0] + '.pdbqt'
        receptor_pt_path = os.path.join(self.tmpdir.name, receptor_pt_path)

        ligand_h_path = self._ligand_path.split('/')[-1].split('.')[0] + 'H.pdb'
        ligand_h_path = os.path.join(self.tmpdir.name, ligand_h_path)

        cmdline_cd = 'cd {}'.format(self.tmpdir.name)

        if not self._receptor_packed: # the packed side chain does not need to add Hs
            self._reduce_receptor(receptor_h_path)
        else:
            receptor_h_path = self._receptor_path

        if not self._ligand_packed:
            self._reduce_ligand(ligand_h_path)
        else:
            ligand_h_path = self._ligand_path

        prep_rec_list = [self.prepare_receptor_bin, '-r', os.path.basename(receptor_h_path)]

        prep_lig_list = [self.prepare_ligand_bin, '-l', os.path.basename(ligand_h_path)]

        cmdline = cmdline_cd + ' && ' + ' '.join(prep_rec_list)
        os.system(cmdline)

        cmdline = cmdline_cd + ' && ' + ' '.join(prep_lig_list)
        os.system(cmdline)

        if auto_box:
            agfr_list = [self.agfr_bin, 
                        '-r',  receptor_h_path.split('.')[0] + '.pdbqt', 
                        '-b', 'user {} {}'.format(' '.join(str(x) for x in docking_center), 
                                                ' '.join(str(x) for x in box_size)),
                        '-l',  ligand_h_path.split('.')[0] + '.pdbqt',
                        '-asv', '1.1',
                        '-o', 'prepared']
        else:
            agfr_list = [self.agfr_bin, 
            '-r',  receptor_h_path.split('.')[0] + '.pdbqt', 
            '-l',  ligand_h_path.split('.')[0] + '.pdbqt',
            '-asv', '1.1',
            '-o', 'prepared']

        cmdline = cmdline_cd + ' && ' + ' '.join(agfr_list)
        os.system(cmdline)

        # Check if directory is empty
        if not os.listdir(self.tmpdir.name):
            print(f"Temporary directory {self.tmpdir.name} is empty - skipping docking")
            return None
            
        # Check for prepared file
        prepared_files = [name for name in os.listdir(self.tmpdir.name) if name[:8] == 'prepared' and name.endswith('.trg')]
        if not prepared_files:
            print(f"No prepared .trg file found in {self.tmpdir.name} - skipping docking") 
            return None
        prepared_file = [name for name in os.listdir(self.tmpdir.name) if name[:8] == 'prepared' and name.endswith('.trg')][0]
        adcp_list = [self.adcp_bin, 
                     '-t', prepared_file,
                     '-s', lig_seq,
                     '-N', str(n_search),
                     '-n', str(n_steps),
                     '-o', 'rl_redocking',
                     '-ref', os.path.basename(self._ligand_path),
                     '>', os.path.join(self.tmpdir.name, 'terminal_record.txt')]
        
        cmdline = cmdline_cd + ' && ' + ' '.join(adcp_list)
        
        os.system(cmdline)

        return self._dump_complex_pdb(save_name, n_save) # click 进去
    
def print_residue_numbers(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", pdb_file)
    
    for model in structure:
        for chain in model:
            print(f"Chain {chain.id} residues:")
            residue_nums = [res.id[1] for res in chain]
            print(residue_nums)


def analyze_hydrogen_bonds(pdb_file, target_residues=[59, 98, 112], max_distance=3.5):
    """Analyze hydrogen bonds between peptide and specific protein residues."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdb_file)
    
    donors = {'N', 'NE', 'NH1', 'NH2', 'ND1', 'NE2', 'NZ', 'OG', 'OG1', 'OH'}
    acceptors = {'O', 'OD1', 'OD2', 'OE1', 'OE2', 'OH'}
    
    h_bonds = []
    
    peptide_chain = structure[0]['X']
    protein_chain = structure[0]['C']  # Adjust chain ID if needed
    
    # Filter protein atoms to only include target residues
    protein_atoms = [atom for residue in protein_chain 
                    if residue.id[1] in target_residues 
                    for atom in residue]
    peptide_atoms = [atom for residue in peptide_chain for atom in residue]
    
    # Check both directions with target residues
    for donor_chain, acceptor_chain, chain_type in [(peptide_atoms, protein_atoms, 'peptide->protein'), 
                                                   (protein_atoms, peptide_atoms, 'protein->peptide')]:
        for donor in donor_chain:
            if donor.name in donors:
                for acceptor in acceptor_chain:
                    if acceptor.name in acceptors:
                        distance = donor - acceptor
                        if distance <= max_distance:
                            donor_res = donor.get_parent()
                            acceptor_res = acceptor.get_parent()
                            
                            h_bonds.append({
                                'type': chain_type,
                                'donor_chain': donor_res.get_parent().id,
                                'donor_residue': f"{donor_res.resname}{donor_res.id[1]}",
                                'donor_atom': donor.name,
                                'acceptor_chain': acceptor_res.get_parent().id,
                                'acceptor_residue': f"{acceptor_res.resname}{acceptor_res.id[1]}",
                                'acceptor_atom': acceptor.name,
                                'distance': round(distance, 2)
                            })
    
    return h_bonds

def print_hbond_analysis(pdb_file, target_residues=[59, 98, 112]):
    h_bonds = analyze_hydrogen_bonds(pdb_file, target_residues)
    
    # Create output string
    output = [f"\nAnalyzing H-bonds for {os.path.basename(pdb_file)} with residues {target_residues}"]
    output.append(f"Found {len(h_bonds)} potential hydrogen bonds:")
    for bond in h_bonds:
        output.append(f"{bond['type']}: {bond['donor_residue']}:{bond['donor_atom']} -> "
                     f"{bond['acceptor_residue']}:{bond['acceptor_atom']} ({bond['distance']}Å)")
    
    # Save to HBOND file in the same directory as the PDB
    hbond_file = os.path.join(os.path.dirname(pdb_file), "HBOND.txt")
    with open(hbond_file, 'a') as f:
        f.write('\n'.join(output) + '\n\n')
    
    # Also print to console
    print('\n'.join(output))
    
if __name__ == '__main__':
    current_dir = os.getcwd()
   
    # optimize_nanomed_ppflow_233k2810/0002_7jjc_2024_10_28__13_31_13/
    protein_name = '1p32'
    optimized_dir = 'codesign_ppflow_233k250218-1M'
    # optimized_peptides = ['0027', '0021', '0033', '0016', '0043', '0020', '0018']
    # optimized_peptides = ['0007', '0045','0046','reference']
    # optimized_peptides = ['reference']
    # optimized_peptides = [f'{i:04d}' for i in range(100)] + ['reference']
    optimized_peptides = ['reference']
    
    protein_file = os.path.join(current_dir, 'dataset', 'nanomed', protein_name, 'receptor.pdb')
    
    save_dir = os.path.join(current_dir, 'results-nanomed','docking', optimized_dir, protein_name)
    os.makedirs(save_dir, exist_ok=True)
    
        
    # Use it before docking
    # print_residue_numbers(protein_file)
    # optimized_peptides = []
    # for peptide in [f'{i:04d}' for i in range(00, 20)]:
    #     bb3_file = os.path.join(current_dir, 'results-nanomed', 'ppflow', 'codesign_nanomed_ppflow_233k/0000_1kex_2024_11_07__11_50_19/', f'{peptide}_bb3.pdb')
    #     if os.path.exists(bb3_file):
    #         optimized_peptides.append(peptide)
            
    for optimized_peptide in optimized_peptides:
        ligand_file = os.path.join(current_dir, 'results-nanomed', 'ppflow', 'codesign_nanomed_ppflow_233k20250217/0002_1p32_2025_02_17__12_07_21', f'{optimized_peptide}.pdb')
        
        
        dock = ADCPDock(save_dir=save_dir)
        dock.set_ligand(ligand_file)
        dock.set_receptor(protein_file) 
        dock.side_chain_packing('ligand')
        np.sum(dock.dock(save_name=f'{protein_name}_{optimized_peptide}_docked_from_{optimized_dir}', n_save=1, auto_box=True))
        
            # Add H-bond analysis
        # docked_pdb = os.path.join(save_dir, f"{protein_name}_{optimized_peptide}_docked_from_{optimized_dir}_1.pdb")
        # print_hbond_analysis(docked_pdb, target_residues=[353, 346, 349])
        