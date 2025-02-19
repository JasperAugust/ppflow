def parse_pdb_file(pdb_file):
   """
   Parses SEQRES lines from a PDB file and returns a dictionary
   mapping chain IDs to their peptide sequence.
   """
   # Try parsing using SEQRES records first
   chain_sequences = {}
   with open(pdb_file, 'r') as f:
       for line in f:
           if line.startswith("SEQRES"):
               tokens = line.split()
               if len(tokens) < 5:
                   continue
               chain_id = tokens[2]
               residues = tokens[4:]
               if chain_id not in chain_sequences:
                   chain_sequences[chain_id] = []
               chain_sequences[chain_id].extend(residues)
   if chain_sequences:  # if SEQRES found, convert residues
       chain_oneletter = {}
       for chain, residues in chain_sequences.items():
           seq = ""
           for res in residues:
               seq += three_to_one.get(res.upper(), "X")
           chain_oneletter[chain] = seq
       return chain_oneletter
   # Fallback: parse ATOM records if no SEQRES found
   chain_residues = {}
   seen = {}
   with open(pdb_file, 'r') as f:
       for line in f:
           if line.startswith("ATOM"):
               chain_id = line[21]
               res_name = line[17:20].strip()
               res_seq = line[22:26].strip()
               key = (chain_id, res_seq)
               if key not in seen:
                   seen[key] = True
                   if chain_id not in chain_residues:
                       chain_residues[chain_id] = []
                   chain_residues[chain_id].append(res_name)
   chain_oneletter = {}
   for chain, residues in chain_residues.items():
       seq = ""
       for res in residues:
           seq += three_to_one.get(res.upper(), "X")
       chain_oneletter[chain] = seq
   return chain_oneletter
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

def __main__():
    parse_pdb_file("peptides.fasta")

if __name__ == "__main__":
    __main__()