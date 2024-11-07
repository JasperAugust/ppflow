import os
import glob
from collections import defaultdict
from Bio import PDB
from Bio.PDB import PDBParser
import pandas as pd

def get_peptide_sequence(pdb_file):
    """Extract peptide sequence from PDB file."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('peptide', pdb_file)
        model = structure[0]
        chain = next(iter(model))
        
        aa_codes = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        
        sequence = ''
        for residue in chain:
            if residue.get_resname() in aa_codes:
                sequence += aa_codes[residue.get_resname()]
        return sequence
    except Exception as e:
        print(f"Error parsing sequence from {pdb_file}: {e}")
        return None

def is_cendr(sequence):
    """Check if sequence ends with R (CendR motif)"""
    return sequence.endswith('R') if sequence else False

def parse_energy_file(file_path):
    """Parse the energy record file and return the best affinity."""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.strip().startswith('1'):
                    return float(line.split()[1]), int(line.split()[4])
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return None, None

def analyze_docking_results(folder_path):
    """Analyze all docking results in the given folder."""
    results = {}
    
    pattern = os.path.join(folder_path, '*_energy_record.txt')
    energy_files = glob.glob(pattern)
    
    for file_path in energy_files:
        peptide_name = os.path.basename(file_path).split('_docked')[0]
        pdb_file = os.path.join(folder_path, f"{peptide_name}_docked_from_optimize_ppflow_233k2810_1.pdb")
        
        sequence = get_peptide_sequence(pdb_file) if os.path.exists(pdb_file) else None
        affinity, cluster_size = parse_energy_file(file_path)
        
        if affinity is not None:
            results[peptide_name] = {
                'affinity': affinity,
                'cluster_size': cluster_size,
                'sequence': sequence,
                'cendr': is_cendr(sequence)
            }
    
    return results

def print_and_save_results(results, output_file):
    """Print results and save to CSV file."""
    print("\nDocking Results:")
    print(f"{'Peptide':<20} {'Sequence':<15} {'Affinity (kcal/mol)':<20} {'Cluster Size':<12} {'CendR':<8}")
    print("-" * 75)
    
    # Sort by affinity
    sorted_results = dict(sorted(results.items(), key=lambda x: x[1]['affinity']))
    
    # Prepare data for DataFrame
    df_data = []
    for peptide, data in sorted_results.items():
        print(f"{peptide:<20} {data['sequence']:<15} {data['affinity']:<20.1f} {data['cluster_size']:<12} {data['cendr']:<8}")
        df_data.append({
            'Peptide': peptide,
            'Sequence': data['sequence'],
            'Affinity': data['affinity'],
            'Cluster_Size': data['cluster_size'],
            'CendR': data['cendr']
        })
    
    # Save to CSV
    df = pd.DataFrame(df_data)
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")

if __name__ == "__main__":
    folder_path = "/gpfs/helios/home/tootsi/homing/ppflow/results-nanomed/docking/optimize_ppflow_233k2810/7jjc"
    output_file = os.path.join(folder_path, "docking_analysis.csv")
    results = analyze_docking_results(folder_path)
    print_and_save_results(results, output_file)