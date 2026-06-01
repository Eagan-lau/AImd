"""
Author: ygliu
Description: This script automates the molecular file preparation process for both protein and ligand types. It incorporates a timeout mechanism to prevent program hangs and allows users to specify additional options for customization. Furthermore, it includes resource management to ensure tasks are executed only when sufficient resources are available and records any commands that exceeded the timeout threshold, facilitating troubleshooting. 
"""

import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import shlex

class ResourceManager:
    def __init__(self, max_attempts=10000, sleep_duration=60):
        self.max_attempts = max_attempts
        self.sleep_duration = sleep_duration

    def has_resources(self):
        try:
            with ThreadPoolExecutor(max_workers=1) as executor:
                executor.submit(print, "Checking resources...")
            return True
        except RuntimeError:
            return False

    def wait_for_resources(self):
        attempts = 0
        while attempts < self.max_attempts:
            if self.has_resources():
                return True
            print("Not enough resources, waiting and retrying...")
            time.sleep(self.sleep_duration)
            attempts += 1
        raise RuntimeError("Max attempts reached. Not enough resources.")

class MoleculePreparation:
    def __init__(self, num_threads=None):
        self.resource_manager = ResourceManager()
        self.num_threads = num_threads
        self.timed_out_commands = []

    def generate_file_list(self, base_path, partition_number):
        target_folder = os.path.join(base_path, f'partition_{partition_number}')
        if not os.path.exists(target_folder):
            print(f"Folder {target_folder} does not exist.")
            return None
        return [f for f in os.listdir(target_folder) if os.path.isfile(os.path.join(target_folder, f))]

    def execute_single_command(self, filename, molecule_type, partition_number, timeout, base_input_path, out_dir, extra_options_ligand, extra_options_receptor, existing_files):
        file_root = os.path.splitext(filename)[0]
        output_filename = f"{file_root}.pdbqt"
        
        # Check if the output file already exists
        if file_root in existing_files:
            print(f"File {output_filename} already exists in {out_dir}, skipping...")
            return
        
        extra_options_receptor = shlex.split(extra_options_receptor)
        extra_options_ligand = shlex.split(extra_options_ligand)
        command1 = []
        command2 = []
        if molecule_type == 'protein':
            command1 = [
                'timeout', str(timeout),
                'prepare_receptor',
                '-r', f'{base_input_path}/partition_{partition_number}/{filename}',
                '-o', f'{out_dir}/{output_filename}',
                '-A', 'hydrogens',
                *extra_options_receptor
            ]
        elif molecule_type == 'ligand':
            if os.path.splitext(filename)[1] not in ['.sdf', '.mol2', '.pdb']:
                print("Only '.sdf', '.mol2', or '.pdb' are supported. Please check the input file type.")
                return
            if os.path.splitext(filename)[1] == '.sdf':
                command1 = [
                    'timeout', str(timeout),
                    'mk_prepare_ligand.py',
                    '-i', f'{base_input_path}/partition_{partition_number}/{filename}',
                    '-o', f'{out_dir}/{output_filename}',
                    *extra_options_ligand
                ]
            elif os.path.splitext(filename)[1] in ['.mol2', '.pdb']:
                command2 = [
                    'timeout', str(timeout),
                    'prepare_ligand',
                    '-l', f'{base_input_path}/partition_{partition_number}/{filename}',
                    '-o', f'{out_dir}/{output_filename}',
                    *extra_options_ligand
                ]
        #print(' '.join(command))
        #subprocess.run(command)
        if command1:
            print(' '.join(command1))
            try:
                subprocess.run(command1, timeout=timeout)
            except subprocess.TimeoutExpired:
                print(f"Command timed out: {' '.join(command1)}")
                self.timed_out_commands.append(' '.join(command1))
                
        if command2:
            print(' '.join(command2))
            try:
                process = subprocess.Popen(command2, cwd=f'{base_input_path}/partition_{partition_number}', stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
            except Exception as e:
                print(f"An error occurred while running the command: {str(e)}")

    def prepare_molecules(self, molecule_type, partition_number, extra_options_ligand, extra_options_receptor, timeout):
        current_dir = os.getcwd()
        base_input_path = f"{current_dir}/data_input/data_{molecule_type}"
        base_output_path = f"{current_dir}/data_output/1_PDBQT/{molecule_type}"

        out_dir = os.path.join(base_output_path, f'file_{partition_number}')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        file_list = self.generate_file_list(base_input_path, partition_number)
        if file_list is None:
            return

        # Read existing output files
        existing_files = [os.path.splitext(f)[0] for f in os.listdir(out_dir) if f.endswith('.pdbqt')]
        
        file_list = self.generate_file_list(base_input_path, partition_number)
        if file_list is None:
            return
            
        self.resource_manager.wait_for_resources()
        with ThreadPoolExecutor(max_workers=self.num_threads) as executor:
            futures = [executor.submit(self.execute_single_command, filename, molecule_type, partition_number, timeout, base_input_path, out_dir, extra_options_ligand, extra_options_receptor, existing_files) for filename in file_list]
            for future in as_completed(futures):
                if future.exception() is not None:
                    print(f"Error in task: {str(future.exception())}")
        
        if self.timed_out_commands:
            with open(f"timeout_{partition_number}.txt", 'w') as f:
                for cmd in self.timed_out_commands:
                    f.write(cmd + '\n')        

if __name__ == '__main__':
    try:
        available_threads = os.cpu_count()
        print(f"Available CPU threads on this machine: {available_threads}")
        
        parser = argparse.ArgumentParser(description="Generate PDBQT format based on molecule type", add_help=False)
        
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
           
        required.add_argument("-mol", required=True, help="Type of molecule: protein or ligand")
        required.add_argument("-num", type=str, required=True, help="File number to specify folder")
        optional.add_argument("--extra_options_ligand", default="", help="Extra options for ligand command (molecule_type=ligand),the parameter content should be enclosed in single quotation marks ('')")
        optional.add_argument("--extra_options_receptor", default="", help="Extra options for receptor command (molecule_type=protein),the parameter content should be enclosed in single quotation marks ('')")
        optional.add_argument("--timeout", type=int, default=3600, help="Time limit in seconds for each command, default is 3600s")

        default_threads = min(32, available_threads)
        optional.add_argument("--threads", type=int, default=default_threads, help=f"Number of threads for parallel computation. Default is {default_threads}")
        
        optional.add_argument("-h", "--help", action="help", help="Display this help message")
        
        args = parser.parse_args()
        print(f"Number of threads used for this task: {args.threads}")
        
        preparation = MoleculePreparation(num_threads=args.threads)
        preparation.prepare_molecules(args.mol, args.num, args.extra_options_ligand, args.extra_options_receptor, args.timeout)
        
        print("The program has successfully concluded its execution.")
    
    except Exception as e:
        print(f"The program encountered an error and ceased execution: {str(e)}")
