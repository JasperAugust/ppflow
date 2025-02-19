import os
import glob
from collections import defaultdict
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
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

def calculate_peptide_properties(sequence):
    """Calculate various peptide properties."""
    if not sequence:
        return {
            'charge': None,
            'is_lipophilic_c_terminus': None,
            'isoelectric_point': None,
            'hydrophobicity': None
        }
    analysis = ProteinAnalysis(sequence)
    charge = analysis.charge_at_pH(7.4)  # Physiological pH
    is_lipophilic_c_terminus = sequence[-1] in ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']  # Example lipophilic residues
    isoelectric_point = analysis.isoelectric_point()
    hydrophobicity = analysis.gravy()  # GRAVY (Grand Average of Hydropathicity) score
    return {
        'charge': charge,
        'is_lipophilic_c_terminus': is_lipophilic_c_terminus,
        'isoelectric_point': isoelectric_point,
        'hydrophobicity': hydrophobicity
    }

def analyze_docking_results(folder_path):
    """Analyze all docking results in the given folder."""
    results = {}
    
    pattern = os.path.join(folder_path, '*_energy_record.txt')
    energy_files = glob.glob(pattern)
    
    for file_path in energy_files:
        peptide_name = os.path.basename(file_path).split('_energy_record.txt')[0]
        pdb_file = os.path.join(folder_path, f"{peptide_name}_1.pdb")
        try:
            sequence = get_peptide_sequence(pdb_file) if os.path.exists(pdb_file) else None
            print(file_path, sequence)
            affinity, cluster_size = parse_energy_file(file_path)
            properties = calculate_peptide_properties(sequence)
            
            if affinity is not None:
                results[peptide_name] = {
                    'affinity': affinity,
                    'cluster_size': cluster_size,
                    'sequence': sequence,
                    'cendr': is_cendr(sequence),
                    'charge': properties['charge'],
                    'is_lipophilic_c_terminus': properties['is_lipophilic_c_terminus'],
                    'isoelectric_point': properties['isoelectric_point'],
                        'hydrophobicity': properties['hydrophobicity']
                    }
        except Exception as e:
            print(f"Error parsing {file_path}: {e}")
    
    return results

def print_and_save_results(results, output_file):
    """Print results and save to CSV file with metadata."""
    print("\nDocking Results:")
    header = f"{'Peptide':<20} {'Sequence':<15} {'Affinity (kcal/mol)':<20} {'CendR':<8} {'Charge':<8} {'Lipophilic C-term':<18} {'pI':<6} {'Hydrophobicity (GRAVY)':<20}"
    print(header)
    print("-" * len(header))
    
    # Sort by affinity
    sorted_results = dict(sorted(results.items(), key=lambda x: x[1]['affinity']))
    
    # Prepare data for DataFrame
    df_data = []
    for peptide, data in sorted_results.items():
        print(f"{peptide:<20} {data['sequence']:<15} {data['affinity']:<20.1f} {data['cendr']:<8} {data['charge']:<8.2f} {data['is_lipophilic_c_terminus']:<18} {data['isoelectric_point']:<6.2f} {data['hydrophobicity']:<20.2f}")
        df_data.append({
            'Peptide': peptide,
            'Sequence': data['sequence'],
            'Affinity (kcal/mol)': data['affinity'],
            'Cluster_Size': data['cluster_size'],
            'CendR motif?': data['cendr'],
            'Charge at pH 7.4': data['charge'],
            'Lipophilic C-terminus': data['is_lipophilic_c_terminus'],
            'Isoelectric_Point': data['isoelectric_point'],
            'Hydrophobicity (GRAVY)': data['hydrophobicity']
        })
    
    # Metadata descriptions
    metadata = [
        "# Docking Analysis Results",
        "# Columns Description:",
        "# Peptide: Identifier of the peptide.",
        "# Sequence: Amino acid sequence.",
        "# Affinity (kcal/mol): Docking binding affinity (more negative means stronger binding).",
        "# Cluster_Size: Size of the cluster from docking results.",
        "# CendR motif?: Boolean indicating if the peptide has a CendR motif (ends with R).",
        "# Charge: Net charge of the peptide at physiological pH (7.4).",
        "# Lipophilic C-terminus: Boolean indicating if the C-terminus is lipophilic.",
        "# Isoelectric_Point: pH at which the peptide has no net charge.",
        "# Hydrophobicity (GRAVY): Grand Average of Hydropathicity score."
    ]
    
    # Save to CSV with metadata as comments
    with open(output_file, 'w') as f:
        for line in metadata:
            f.write(line + '\n')
        df = pd.DataFrame(df_data)
        df.to_csv(f, index=False)
    
    print(f"\nResults saved to: {output_file}")

if __name__ == "__main__":
    folder_path = "/gpfs/helios/home/tootsi/homing/ppflow/results-nanomed/docking/codesign_ppflow_233k250218-1M/1p32"
    output_file = os.path.join(folder_path, "docking_analysis.csv")
    results = analyze_docking_results(folder_path)
    print_and_save_results(results, output_file)
    