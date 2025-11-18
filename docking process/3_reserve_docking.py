"""
Author: ygliu
Description: This script automates molecular docking simulations with AutoDock Vina. It conducts docking configurations in a specified directory (identified by 'num') and allows optional log file storage ('log_dir'). The script leverages parallel processing, enabling efficient multi-threaded execution. Users can customize parameters like exhaustiveness and additional options. A timeout mechanism prevents potential stalls, and timed-out commands are logged for analysis.
20240527
"""

import os
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def find_pdbqt_files(file_dir, extension='.txt'):
    file_list = []
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if file.endswith(extension):
                file_name_without_extension = os.path.splitext(file)[0]
                file_list.append(file_name_without_extension)
    return file_list

def generate_existing_files(dire, extension='.out'):
    existing_files = []
    for root, dirs, files in os.walk(dire):
        for file in files:
            if file.endswith(extension):
                name = os.path.splitext(file)[0]
                parts = name.split('_cavity_')
                if len(parts) == 2:
                    num = parts[1]
                    base_name = parts[0].replace('_', '@', 1)
                    output_name = f"{base_name}_c{num}_out.pdbqt"
                    output_file_path = os.path.join(root, output_name)
                    if os.path.exists(output_file_path):
                        existing_files.append(name)
    return existing_files

class RunCmd:
    @staticmethod
    def cmd_run(cmd):
        subprocess.call(cmd, shell=True)

timed_out_commands = []
def run_with_timeout(command, timeout):
    try:
        process = subprocess.Popen(command, shell=True)
        process.communicate(timeout=timeout)
        return process.returncode
    except subprocess.TimeoutExpired:
        print(f"Timeout expired for command: {command}")
        timed_out_commands.append(command)
        process.terminate()
        return None

def vina_docking(config_file, output_file, exhaustiveness, additional_params, timeout, cpu, existing_files):
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    file_name_without_extension = os.path.splitext(os.path.basename(config_file))[0]
    if existing_files is not None:
        if file_name_without_extension in existing_files:
            print(f"Skipping {file_name_without_extension} as it already exists.")
            return None
    command = f"timeout {timeout} vina --config {config_file} --exhaustiveness {exhaustiveness} --cpu {cpu} {additional_params} > {output_file}"
    print(command)
    return run_with_timeout(command, timeout)

def main(log_dir=None, num=None, threads=None, exhaustiveness=8, additional_params="", timeout=None, cpu=None, overwrite=False):
    available_cpus = os.cpu_count()
    print(f"Available CPU threads on this machine: {available_cpus}")
    
    if cpu is None or cpu > available_cpus:
        cpu = min(exhaustiveness, available_cpus)
    print(f"CPUs used per vina process: {cpu}")
    
    vina_threads = max(1, int(threads // cpu))
    print(f"Threads used for docking tasks: {vina_threads}") 
    
    current_path = os.getcwd()
    dire_conf = f'{current_path}/data_output/3_pocket_parameters/file_{num}'
    dire = f'{current_path}/data_output/4_docking_results/file_{num}'

    files_to_run = find_pdbqt_files(dire_conf, extension='.txt')
    
    existing_files = []
    if not overwrite:
        existing_files = generate_existing_files(dire, extension='.out')
        #print(existing_files)
    else:
        existing_files = None
    
    with ThreadPoolExecutor(max_workers=vina_threads) as executor:
        futures = []
        for file in files_to_run:
            config_file = f"{dire_conf}/{file}.txt"
            if log_dir:
                output_file = f"{log_dir}/{file}.out"
            else:
                output_file = f"{dire}/{file}.out"
            futures.append(executor.submit(vina_docking, config_file, output_file, exhaustiveness, additional_params, timeout, cpu, existing_files))

        for future in as_completed(futures):
            if future.exception():
                print(f"Error: {future.exception()}")
                
        if timed_out_commands:
            with open(f"timeout_{num}.txt", "a") as f:
                for cmd in timed_out_commands:
                    f.write(cmd + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automating molecular docking simulations", add_help=False)
    
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('-num', type=int, required=True, help="Number for folder identification.")
    optional.add_argument('--log_dir', type=str, help="Directory where the log files will be stored.")
    
    available_threads = os.cpu_count()
    default_threads = min(64, available_threads)
    optional.add_argument('--threads', type=int, default=default_threads, help=f"Number of threads for computation. Default is {default_threads}")

    optional.add_argument("--exhaustiveness", type=int, default=8, help="Exhaustiveness, default is 8")
    optional.add_argument("--cpu", type=int, default=None, help="Number of CPUs for each vina command to use, default is dynamically determined") 
    optional.add_argument("--additional_params", default="", help="Additional parameters, the parameter content should be enclosed in single quotation marks ('')")
    optional.add_argument("--timeout", type=int, default=3600, help="Timeout for the command in seconds, default is 36000s")
    optional.add_argument('--overwrite', action='store_true', help="If specified, the existing output files will be overwritten ")
    
    optional.add_argument("-h", "--help", action="help", help="Display this help message")
    
    args = parser.parse_args()
    
    main(log_dir=args.log_dir, num=args.num, threads=args.threads, exhaustiveness=args.exhaustiveness, additional_params=args.additional_params, timeout=args.timeout, cpu=args.cpu, overwrite=args.overwrite)
    
    print("Molecular docking simulations have been completed")
