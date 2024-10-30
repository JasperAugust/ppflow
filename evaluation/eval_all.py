import argparse
import os 
import numpy as np
import threading
import concurrent.futures
import subprocess

def run_command(command):
    process = subprocess.Popen(command, shell=True)
    process.wait()

def run_commands(commands, block=4):
    with concurrent.futures.ThreadPoolExecutor(max_workers=block) as executor:
        for i in range(0, len(commands), block):
            chunk = commands[i:i+block]
            futures = [executor.submit(run_command, cmd) for cmd in chunk]
            concurrent.futures.wait(futures)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta_gen_dir', type=str, default='/gpfs/helios/home/tootsi/homing/ppflow/results-nanomed/docking/optimize_ppflow_233k2710')
    parser.add_argument('--benchmark_dir', type=str, default='/gpfs/helios/home/tootsi/homing/ppflow/dataset/nanomed')
    parser.add_argument('--mode', type=str, choices=['basic_prop', 'bind_aff'], default='basic_prop')
    parser.add_argument('--parallel', action='store_true', default=False, help='Enable parallel execution')
    parser.add_argument('--no-parallel', action='store_false', dest='parallel', help='Disable parallel execution')
    parser.add_argument('--block', type=int, default=4)
    parser.add_argument('--num_samples', type=int, default=-1, help='Number of samples to evaluate. -1 means all samples.')

    args = parser.parse_args()

    gen_dir_names = sorted(os.listdir(args.meta_gen_dir))


    cmds = []
    for gen_dir_name in gen_dir_names:
        gen_dir = os.path.join(args.meta_gen_dir, gen_dir_name)
        try:
            pdb_id = gen_dir_name.split('_')[1]
        except IndexError:
            print(f"Skipping invalid directory name format: {gen_dir_name}")
            continue
        ref_dir = os.path.join(args.benchmark_dir, pdb_id)

        if args.mode == 'basic_prop':
            cmd = f'python eval_struct_seq.py --gen_dir "{gen_dir}" --ref_dir "{ref_dir}" --save_path "{gen_dir}" --num_samples {args.num_samples}'
            cmds.append(cmd)
        
        elif args.mode == 'bind_aff':
            cmd = f'python eval_bind.py --gen_dir "{gen_dir}" --ref_dir "{ref_dir}" --save_path "{gen_dir}" --num_samples {args.num_samples}'
            cmds.append(cmd)

    if args.parallel:
        run_commands(cmds, args.block)
    else:
        for cmd in cmds:
            os.system(cmd)