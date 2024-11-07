import os
import sys
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



pdb_parser = PDBParser(QUIET=True)


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
    def set_receptor(self, pdb_path):
        pass

    @abc.abstractmethod
    def set_ligand(self, pdb_path):
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


import os
import subprocess
import tempfile
import logging

class ADCPDock(DockingEngine):
    def __init__(
        self, 
        save_dir,
        adcp_suite_path='/gpfs/helios/home/tootsi/homing/ppflow/bin/ADFRsuite_x86_64Linux_1.0',
    ):
        super().__init__()
        self.adcp_suite_path = adcp_suite_path
        self.bin_path = os.path.join(adcp_suite_path, 'bin')
        self.lib_path = os.path.join(adcp_suite_path, 'lib')
        self.python2_bin = os.path.join(self.bin_path, 'python')
        self.tmpdir = tempfile.TemporaryDirectory()
        self.save_dir = save_dir
        os.makedirs(self.save_dir, exist_ok=True)

        self._has_receptor = False
        self._has_ligand = False
        self._receptor_packed = False
        self._ligand_packed = False

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
                        os.path.join(self.save_dir, "{}_{}_energy_record.txt".format(save_name, n_save)))
        
        energies = []
        with open(os.path.join(self.tmpdir.name, "terminal_record.txt"), 'r') as f:
            data = f.read().splitlines()
            
            energs = []

            for i, line in enumerate(data):

                if line and line[0] == "-":
                    energs = data[i+1 : -1]

            for e in energs:

                energies.append(float(e.split()[1]))

        if len(energies) == 0:
            print('Failed!')
            return np.array([np.nan])

        return np.array(energies)

        

    def load_pep_seq_and_struc(self):
        structure_ligand = pdb_parser.get_structure(self._ligand_path, self._ligand_path)[0]

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

    def auto_dock_box(self, lig_struc):
        docking_center = lig_struc.mean(0)
        box_size = np.abs(lig_struc - docking_center[None,:]).max(0) * 1.5
        return docking_center, box_size

    def _run_adcp_command(self, command, *args):
        full_command = [self.python2_bin, os.path.join(self.bin_path, command)] + list(args)
        env = os.environ.copy()
        env['LD_LIBRARY_PATH'] = f"{self.lib_path}:{env.get('LD_LIBRARY_PATH', '')}"
        result = subprocess.run(full_command, env=env, check=True, capture_output=True, text=True)
        return result.stdout

    def dock(self, save_name, n_save, n_search=100, n_steps=100000, auto_box=False):
        logger = logging.getLogger(__name__)

        if not (self._has_receptor and self._has_ligand):
            raise ValueError('Missing receptor or peptide.')

        if auto_box:
            lig_seq, lig_struc = self.load_pep_seq_and_struc()
            docking_center, box_size = self.auto_dock_box(lig_struc)

        receptor_h_path = os.path.join(self.tmpdir.name, f"{os.path.splitext(os.path.basename(self._receptor_path))[0]}H.pdb")
        ligand_h_path = os.path.join(self.tmpdir.name, f"{os.path.splitext(os.path.basename(self._ligand_path))[0]}H.pdb")

        if not self._receptor_packed:
            logger.info(f"Reducing receptor: {self._receptor_path}")
            self._run_adcp_command('reduce', '-HIS', self._receptor_path, '>', receptor_h_path)
        else:
            receptor_h_path = self._receptor_path

        if not self._ligand_packed:
            logger.info(f"Reducing ligand: {self._ligand_path}")
            self._run_adcp_command('reduce', '-HIS', self._ligand_path, '>', ligand_h_path)
        else:
            ligand_h_path = self._ligand_path

        logger.info(f"Preparing receptor: {receptor_h_path}")
        self._run_adcp_command('prepare_receptor', '-r', os.path.basename(receptor_h_path))

        logger.info(f"Preparing ligand: {ligand_h_path}")
        self._run_adcp_command('prepare_ligand', '-l', os.path.basename(ligand_h_path), '-A', 'hydrogens')

        if auto_box:
            agfr_args = ['-r', f"{os.path.splitext(receptor_h_path)[0]}.pdbqt",
                         '-b', f"user {' '.join(map(str, docking_center))} {' '.join(map(str, box_size))}",
                         '-l', f"{os.path.splitext(ligand_h_path)[0]}.pdbqt",
                         '-asv', '1.1',
                         '-o', 'prepared']
        else:
            agfr_args = ['-r', f"{os.path.splitext(receptor_h_path)[0]}.pdbqt",
                         '-l', f"{os.path.splitext(ligand_h_path)[0]}.pdbqt",
                         '-asv', '1.1',
                         '-o', 'prepared']

        logger.info(f"Running AGFR: {' '.join(['agfr'] + agfr_args)}")
        self._run_adcp_command('agfr', *agfr_args)

        prepared_file = next((name for name in os.listdir(self.tmpdir.name) if name.startswith('prepared') and name.endswith('.trg')), None)
        if not prepared_file:
            raise FileNotFoundError(f"No prepared .trg files found in {self.tmpdir.name}")

        logger.info(f"Running ADCP docking")
        adcp_args = ['-t', prepared_file,
                     '-s', lig_seq,
                     '-N', str(n_search),
                     '-n', str(n_steps),
                     '-o', 'rl_redocking',
                     '-ref', os.path.basename(self._ligand_path)]
        self._run_adcp_command('adcp', *adcp_args)

        return self._dump_complex_pdb(save_name, n_save)
    
# if __name__ == '__main__':
#     # protein_name = '2qbx'
#     # ligand_variant = '0006'
#     # protein_file = f'/gpfs/helios/home/tootsi/homing/ppflow/dataset/PPDbench/{protein_name}/receptor.pdb'
#     # ligand_file = f'/gpfs/helios/home/tootsi/homing/ppflow/results/ppflow/codesign_ppflow/0000_{protein_name}_2024_08_19__11_28_55/{ligand_variant}_bb3.pdb'
#     # # ligand_file= f'/gpfs/helios/home/tootsi/homing/ppflow/results/ppflow/codesign_ppflow/0000_{protein_name}_2024_08_19__11_28_55/reference.pdb'
    
#     # results_dir = os.path.join('/gpfs/helios/home/tootsi/homing/ppflow/results-jasper/docking', protein_name)
#     # os.makedirs(results_dir, exist_ok=True)
    
    
    
#     save_name = f'{protein_name}-{ligand_variant}'
#     dock = ADCPDock(save_dir=results_dir)
#     dock.set_ligand(ligand_file)
#     dock.set_receptor(protein_file) 
#     dock.side_chain_packing('ligand')
#     energy_sum = np.sum(dock.dock(save_name=save_name, n_save=1, auto_box=True))
    
#     print(f'Energy sum: {energy_sum}')



if __name__ == '__main__':
    import logging
    from collections import defaultdict
    
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    base_dir = '/gpfs/helios/home/tootsi/homing/ppflow'
    optimize_dir = os.path.join(base_dir, 'results-jasper', 'ppflow', 'optimize_ppflow_233k1410')

    # Group peptides by receptor
    receptor_peptides = defaultdict(list)
    for folder in os.listdir(optimize_dir):
        receptor_name = folder.split('_')[1]
        receptor_peptides[receptor_name].append(folder)

    for receptor_name, folders in receptor_peptides.items():
        # Sort folders to get the last one
        last_folder = sorted(folders)[-1]
        results_dir = os.path.join(optimize_dir, last_folder)
        protein_file = os.path.join(base_dir, 'dataset', 'PPDbench', receptor_name, 'receptor.pdb')

        if not os.path.exists(protein_file):
            logger.error(f"Receptor file not found: {protein_file}")
            continue

        # Get the last peptide file (excluding reference.pdb)
        peptide_files = [f for f in os.listdir(results_dir) if f.endswith('.pdb') and f != 'reference.pdb']
        if not peptide_files:
            logger.warning(f"No peptide files found in {results_dir}")
            continue
        
        last_peptide = sorted(peptide_files)[-1]
        ligand_file = os.path.join(results_dir, last_peptide)
        save_name = f'{receptor_name}-{last_peptide.split(".")[0]}'

        try:
            dock = ADCPDock(save_dir=results_dir)
            dock.set_ligand(ligand_file)
            dock.set_receptor(protein_file)
            dock.side_chain_packing('ligand')
            
            logger.info(f"Starting docking for {save_name}")
            energy_sum = np.sum(dock.dock(save_name=save_name, n_save=1, auto_box=True))
            logger.info(f'{save_name} - Energy sum: {energy_sum}')
        except Exception as e:
            logger.error(f'Error processing {save_name}: {str(e)}', exc_info=True)

    logger.info('ADCP docking completed for the last optimized peptide of each receptor.')