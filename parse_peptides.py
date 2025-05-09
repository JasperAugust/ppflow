#!/usr/bin/env python3
import os
import glob

# folder containing .pdb files
pdb_folder = "/gpfs/helios/home/tootsi/homing/ppflow/results-nanomed/ppflow/codesign_nanomed_ppflow_233k20250217/0002_1p32_2025_02_17__12_07_21"
output_fasta = "peptides20250304.fasta"

# mapping from three-letter to one-letter amino acid codes
three_to_one = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # add more mappings if needed (e.g. SEC, PYL, etc.)
}

def parse_pdb_file(pdb_file):
    """
    Parses a PDB file to extract peptide sequences by first checking SEQRES records.
    If none found, falls back to parsing ATOM records. Debug logs are added to validate the process.
    """
    print("DEBUG: Parsing file", pdb_file)

    # Try parsing using SEQRES records first
    chain_sequences = {}
    seqres_count = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("SEQRES"):
                seqres_count += 1
                tokens = line.split()
                if len(tokens) < 5:
                    continue
                chain_id = tokens[2]
                residues = tokens[4:]
                if chain_id not in chain_sequences:
                    chain_sequences[chain_id] = []
                chain_sequences[chain_id].extend(residues)
    print(f"DEBUG: Found {seqres_count} SEQRES records in {pdb_file}")

    if chain_sequences:  # if SEQRES found, convert residues
        for chain, residues in chain_sequences.items():
            print(f"DEBUG: Chain {chain} has {len(residues)} residues from SEQRES in {pdb_file}")
        chain_oneletter = {}
        for chain, residues in chain_sequences.items():
            seq = ""
            for res in residues:
                seq += three_to_one.get(res.upper(), "X")
            chain_oneletter[chain] = seq
        return chain_oneletter

    # Fallback: parse ATOM records if no SEQRES found
    print(f"DEBUG: No SEQRES records found in {pdb_file}. Using ATOM records as fallback.")
    chain_residues = {}
    seen = {}
    atom_count = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                atom_count += 1
                chain_id = line[21]
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                key = (chain_id, res_seq)
                if key not in seen:
                    seen[key] = True
                    if chain_id not in chain_residues:
                        chain_residues[chain_id] = []
                    chain_residues[chain_id].append(res_name)
    print(f"DEBUG: Processed {atom_count} ATOM lines in {pdb_file}. Found {len(chain_residues)} chains.")

    chain_oneletter = {}
    for chain, residues in chain_residues.items():
        print(f"DEBUG: Chain {chain} has {len(residues)} residues from ATOM fallback in {pdb_file}")
        seq = ""
        for res in residues:
            seq += three_to_one.get(res.upper(), "X")
        chain_oneletter[chain] = seq
    return chain_oneletter

def main():
    pdb_files = [f for f in glob.glob(os.path.join(pdb_folder, "*[pP][dD][bB]")) if not f.endswith("_bb3.pdb")]
    if not pdb_files:
        print("No PDB files found in folder:", pdb_folder)
        return

    with open(output_fasta, "w") as outfile:
        for pdb_file in pdb_files:
            print("Processing", pdb_file)
            sequences = parse_pdb_file(pdb_file)
            if not sequences:
                print("No sequences found in", pdb_file)
                continue
            filename = os.path.basename(pdb_file)
            for chain, seq in sequences.items():
                header = f">{filename}_{chain}"
                outfile.write(header + "\n")
                for i in range(0, len(seq), 60):
                    outfile.write(seq[i:i+60] + "\n")
    print(f"FASTA file generated: {output_fasta}")

if __name__ == "__main__":
    main()