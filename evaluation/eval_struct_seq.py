import argparse
from tools.score import rosetta_energy, similarity, check_validity, foldx_energy
import os 
import numpy as np
import logging
import pandas as pd
import shutil

logging.basicConfig(level=logging.INFO) 
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gen_dir', type=str, default='./results/diffpp/codesign_diffpp/0000_2qbx_2024_01_18__18_37_25')
    parser.add_argument('--ref_dir', type=str, default='./PPDbench/2qbx/')
    parser.add_argument('--save_path', type=str, default='./results/diffpp/codesign_diffpp/0000_2qbx_2024_01_18__18_37_25')
    parser.add_argument('--num_samples', type=int, default=-1, help='Number of samples to evaluate. -1 means all samples.')
    args = parser.parse_args()
    
    logger.info('Evaluating the samples in {}'.format(args.gen_dir))
    pdb_names = sorted([f for f in os.listdir(args.gen_dir) if f.endswith('.pdb')])
    if args.num_samples > 0:
        pdb_names = pdb_names[:args.num_samples]
    valid_pdb_names = []

    for pdb_name in pdb_names:
        try:
            valid_bool = check_validity.bond_length_validation(os.path.join(args.gen_dir, pdb_name))
            if valid_bool:
                valid_pdb_names.append(pdb_name)
        except Exception as e:
            logger.error(f"Error validating {pdb_name}: {e}")
    
    metrics = {'dg': [], 'novel': [], 'seq_div': [], 'str_div': []}
    metrics['validity'] = len(valid_pdb_names) / len(pdb_names) if pdb_names else 0
    ref_file = os.path.join(args.ref_dir, 'peptide.pdb')
    protein_file = os.path.join(args.ref_dir, 'receptor.pdb')
    
    try:
        score_calculator = foldx_energy.FoldXGibbsEnergy()
        score_calculator.set_ligand(ref_file)
        score_calculator.set_receptor(protein_file)
        score_calculator.side_chain_packing('ligand')
        dg_ref = score_calculator.cal_interface_energy()
        metrics['dg_ref'] = dg_ref
    except Exception as e:
        logger.error(f"Error calculating reference energy: {e}")
        metrics['dg_ref'] = None
    
    # Initialize variables to track minimum energy
    min_dg = float('inf')
    min_dg_peptide = None

    # Initialize list to collect per-peptide metrics
    peptide_metrics = []

    for pdb_name in valid_pdb_names:
        try:
            score_calculator = foldx_energy.FoldXGibbsEnergy()
            gen_bb_file = os.path.join(args.gen_dir, pdb_name)
            score_calculator.set_ligand(gen_bb_file)
            score_calculator.set_receptor(protein_file)
            score_calculator.side_chain_packing('ligand')
            dg = score_calculator.cal_interface_energy()
            metrics['dg'].append(dg)

            # Check and update minimum energy
            if dg < min_dg:
                min_dg = dg
                min_dg_peptide = pdb_name

            novel_bool = similarity.check_novelty(gen_bb_file, ref_file)
            metrics['novel'].append(novel_bool)

            for pdb_name_vs in valid_pdb_names:
                vs_pdb_file = os.path.join(args.gen_dir, pdb_name_vs)

                seq_sim = similarity.seq_similarity(vs_pdb_file, gen_bb_file)
                metrics['seq_div'].append(1 - seq_sim)

                str_sim = similarity.tm_score(vs_pdb_file, gen_bb_file)
                metrics['str_div'].append(1 - str_sim)
                
            # Append per-peptide metrics
            peptide_metrics.append({
                'pdb_name': pdb_name,
                'dg': dg,
                'novel': novel_bool,
                'seq_div': 1 - seq_sim,
                'str_div': 1 - str_sim
            })

        except Exception as e:
            logger.error(f"Error processing {pdb_name}: {e}")
    
    # Add minimum energy metrics
    metrics['min_dg'] = min_dg
    metrics['min_dg_peptide'] = min_dg_peptide

    # Create a DataFrame from the collected metrics
    df_metrics = pd.DataFrame(peptide_metrics)

    # Define the path for the CSV file
    csv_path = os.path.join(args.save_path, 'peptide_metrics.csv')

    # Save the DataFrame to CSV
    df_metrics.to_csv(csv_path, index=False)
    logger.info(f"Saved peptide metrics to {csv_path}")

    np.save(os.path.join(args.save_path, 'metrics_meta'), metrics)
    logger.info('Saved metrics to {}'.format(os.path.join(args.save_path, 'metrics_meta.npy')))