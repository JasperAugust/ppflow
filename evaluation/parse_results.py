import os
import numpy as np
import pandas as pd

def parse_results(results_dir):
    metrics = {}
    for folder in os.listdir(results_dir):
        folder_path = os.path.join(results_dir, folder)
        if os.path.isdir(folder_path):
            metrics_file = os.path.join(folder_path, 'metrics_meta.npy')
            if os.path.exists(metrics_file):
                data = np.load(metrics_file, allow_pickle=True).item()
                pdb_id = folder.split('_')[1]
                metrics[pdb_id] = {
                    'validity': data['validity'],
                    'dg_ref': data['dg_ref'],
                    'dg_mean': np.mean(data['dg']),
                    'dg_std': np.std(data['dg']),
                    'novelty': np.mean(data['novel']),
                    'seq_div_mean': np.mean(data['seq_div']),
                    'str_div_mean': np.mean(data['str_div'])
                }
    
    df = pd.DataFrame.from_dict(metrics, orient='index')
    return df

if __name__ == '__main__':
    results_dir = '../results-jasper/ppflow/optimize_ppflow_233k2410/'
    df = parse_results(results_dir)
    print(df)
    df.to_csv('../results-jasper/ppflow/optimize_ppflow_233k2410/parsed_results_optimize_ppflow_233k2410.csv')
