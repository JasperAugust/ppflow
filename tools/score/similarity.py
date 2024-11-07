#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import re
import time
from ppflow.datasets.constants import *
import numpy as np 
from Bio.PDB import PDBIO, PDBParser
parser = PDBParser(QUIET=True)

FILE_DIR = os.path.split(__file__)[0]
TMEXEC = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'bin', 'TMscore', 'TMscore')

def load_pep_seq(pdb_path):
    structure_ligand = parser.get_structure(pdb_path, pdb_path)[0]

    for chain in structure_ligand:
        seq = []
        for residue in chain:
            seq.append(resindex_to_ressymb[AA(residue.get_resname())])
    seq = ''.join(seq)
    seq = seq.lower()

    return seq


def tm_score(pdb_path1, pdb_path2):
    p = os.popen(f'{TMEXEC} {pdb_path1} {pdb_path2}')
    text = p.read()
    p.close()
    res = re.search(r'TM-score\s*= ([0-1]\.[0-9]+)', text)
    score = float(res.group(1))
    return score

def seq_similarity(pdb_path1, pdb_path2):
    seq1 = load_pep_seq(pdb_path1)
    seq2 = load_pep_seq(pdb_path2)
    list1 = np.array(list(seq1))
    list2 = np.array(list(seq2))
    similarity = (list1 == list2).sum() / len(list1)
    return similarity

def check_novelty(pdb_path1, pdb_path2, tm_threshold=0.5, seq_threshold=0.5):
    tm = tm_score(pdb_path1, pdb_path2)
    seq_sim = seq_similarity(pdb_path1, pdb_path2)
    if tm < tm_threshold and seq_sim < seq_threshold:
        return True
    return False


if __name__ == '__main__':
    # TODO note(Jasper): Replace this with a specific file!
    import os
    import pandas as pd
    codesign_ppflow_dir = '/gpfs/helios/home/tootsi/homing/ppflow/results/ppflow/codesign_ppflow'
    
    for folder in os.listdir(codesign_ppflow_dir):
        base_dir = os.path.join(codesign_ppflow_dir, folder)
        if not os.path.isdir(base_dir):
            continue
        
        pdb_id = folder.split('_')[1]
        pdb_path2 = f'/gpfs/helios/home/tootsi/homing/ppflow/dataset/PPDbench/{pdb_id}/peptide.pdb'
        
        results = []
        
        for filename in os.listdir(base_dir):
            if filename.endswith('.pdb'):
                pdb_path1 = os.path.join(base_dir, filename)
                tm = tm_score(pdb_path1, pdb_path2)
                seq_sim = seq_similarity(pdb_path1, pdb_path2)
                novelty = check_novelty(pdb_path1, pdb_path2)
                
                results.append({
                    'Filename': filename,
                    'TM Score': tm,
                    'Sequence Similarity': seq_sim,
                    'Novelty': novelty
                })
        
        df = pd.DataFrame(results)
        
        # Save to file
        output_file = os.path.join(base_dir, 'similarity_results.csv')
        df.to_csv(output_file, index=False)
        
        # Print to console
        print(f"\nResults for {folder}:")
        print(df.to_string(index=False))
        print(f"Results saved to: {output_file}")