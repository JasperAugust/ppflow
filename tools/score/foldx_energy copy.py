import tempfile
import shutil
import os
from tools.relax.rosetta_packing import side_chain_packing
from tools.base import merge_protein_ligand

class FoldXGibbsEnergy():

    def __init__(
        self,
        foldx_path = './bin/FoldX/foldx_20241231'
    ):
        super().__init__()
        self.foldx_path = os.path.abspath(foldx_path)
        self.tmpdir = tempfile.TemporaryDirectory()

        self._has_receptor = False
        self._has_ligand = False

        self._receptor_chains = []
        self._ligand_chains = []


    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        self.tmpdir.cleanup()

    def set_receptor(self, pdb_path):
        shutil.copyfile(pdb_path, os.path.join(self.tmpdir.name, 'receptor.pdb'))
        self._has_receptor = True
        self._receptor_path = os.path.join(self.tmpdir.name, 'receptor.pdb')

    def set_ligand(self, pdb_path):
        shutil.copyfile(pdb_path, os.path.join(self.tmpdir.name, 'ligand.pdb'))
        self._has_ligand = True
        self._ligand_path = os.path.join(self.tmpdir.name, 'ligand.pdb')

    def side_chain_packing(self, type):
        out_path = side_chain_packing(pdb_file=os.path.join(self.tmpdir.name, '{}.pdb'.format(type)),
                                      output_file=os.path.join(self.tmpdir.name, '{}_packed.pdb'.format(type)))
        
        if type == 'ligand':
            self._ligand_path = out_path
        elif type == 'receptor':
            self._receptor_path = out_path
        else:
            raise ValueError('No such type')

        return out_path
    
    def merge_protein_ligand(self):
        out_file, self.interface = merge_protein_ligand(self._receptor_path, self._ligand_path, 
                                        out_pdb_file=os.path.join(self.tmpdir.name, 'merge.pdb'))
        return out_file
        
    def cal_interface_energy(self):
        merge_file = self.merge_protein_ligand()
        cmd="cd "+ self.tmpdir.name +"; "
        cmd += self.foldx_path + " --command=Stability" + " --pdb="+ os.path.basename(merge_file)

        os.system(cmd)

        output_file = os.path.join(self.tmpdir.name, 'merge_0_ST.fxout')
        if not os.path.exists(output_file):
            raise FileNotFoundError(f"FoldX output file not found: {output_file}")
        return self.read_output(output_file)

    def read_output(self, file_name):
        with open(file_name, 'r') as f:
            data = f.read().splitlines()
            total_energy = float(data[0].split('\t')[1])
        return total_energy



# Jasper 
# if __name__ == '__main__':
#     import os
#     import pandas as pd
    
#     codesign_ppflow_dir = '/gpfs/helios/home/tootsi/homing/ppflow/results/ppflow/codesign_ppflow'
    
#     for folder in os.listdir(codesign_ppflow_dir):
#         base_dir = os.path.join(codesign_ppflow_dir, folder)
#         if not os.path.isdir(base_dir):
#             continue
        
#         pdb_id = folder.split('_')[1]
#         protein_file = f'/gpfs/helios/home/tootsi/homing/ppflow/dataset/PPDbench/{pdb_id}/receptor.pdb'
        
#         results = []
        
#         for filename in os.listdir(base_dir):
#             if filename.endswith('.pdb'):
#                 gen_bb_file = os.path.join(base_dir, filename)
                
#                 score_calculator = FoldXGibbsEnergy()
#                 score_calculator.set_ligand(gen_bb_file)
#                 score_calculator.set_receptor(protein_file)
#                 score_calculator.side_chain_packing('ligand')
#                 dg = score_calculator.cal_interface_energy()
                
#                 results.append({
#                     'Filename': filename,
#                     'FoldX Energy': dg
#                 })
        
#         df = pd.DataFrame(results)
        
#         # Save to file
#         output_file = os.path.join(base_dir, 'foldx_energy_results.csv')
#         df.to_csv(output_file, index=False)
        
#         # Print to console
#         print(f"\nResults for {folder}:")
#         print(df.to_string(index=False))
#         print(f"Results saved to: {output_file}")


