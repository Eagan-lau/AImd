import os
import argparse
from multiprocessing import Pool

def read_data_from_file(file_path, num_affinities):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            affinities = []
            grid_size = grid_space = exhaustiveness = random_seed = None
            
            for line in lines:
                if line.startswith("Grid size"):
                    parts = line.split(':')[1].strip().split()
                    grid_size = ",".join(parts[i] for i in range(1, len(parts), 2))  # Format as "34,28,30"
                elif line.startswith("Grid space"):
                    grid_space = line.split(':')[1].strip()
                elif line.startswith("Exhaustiveness"):
                    exhaustiveness = line.split(':')[1].strip()
                elif line.startswith("Performing docking"):
                    random_seed = line.split('(')[1].split(')')[0].strip().replace('random seed: ', '')  # Get value inside parentheses

                # Check for affinity lines based on the specified number
                if len(affinities) < num_affinities:
                    affinity_line_index = len(affinities) + 1
                    if line.startswith(f"   {affinity_line_index}"):
                        parts = line.split()
                        affinities.append(float(parts[1]))

            # Fill with 'NA' if there are not enough affinity values
            while len(affinities) < num_affinities:
                affinities.append('NA')
                
            return affinities, grid_size, grid_space, exhaustiveness, random_seed
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, None, None, None, None

def process_file(file_path, num_affinities):
    file_name = os.path.basename(file_path)
    result = read_data_from_file(file_path, num_affinities)
    if result[0] is not None:  # Check if we successfully read data
        return (file_name, *result)
    return (file_name, [], None, None, None, None)

def find_files_in_dir(file_dir, extension='.out'):
    file_list = []
    for root, dirs, files in os.walk(file_dir):
        for file in files:
            if file.endswith(extension):
                file_path = os.path.join(root, file)
                file_list.append(file_path)
    return file_list

def main(folder_range=None, num_affinities=1, num_threads=None):
    results = {}
    current_path = os.getcwd()
    parent_directory = f'{current_path}/data_output/4_docking_results'
    output_directory = f'{current_path}/data_output/5_statistics'
    
    os.makedirs(output_directory, exist_ok=True)  # Create output directory if it doesn't exist
    
    if folder_range:
        try:
            start_num, end_num = map(int, folder_range.split('-'))
            folder_range = range(start_num, end_num + 1)
        except ValueError:
            print("Invalid folder range format. Please use 'start-end'.")
            return
    else:
        folder_range = [int(folder.split('_')[-1]) for folder in os.listdir(parent_directory) if folder.startswith('file_')]
    
    file_paths = []
    for num in folder_range:
        directory = f'{parent_directory}/file_{num}'
        file_paths.extend(find_files_in_dir(directory))
    
    # Automatically detect the number of available threads
    available_threads = os.cpu_count()
    if available_threads > 32:
        available_threads = 32
    if num_threads:
        available_threads = min(available_threads, int(num_threads))
    
    print(f"Using {available_threads} threads for processing.")
    
    with Pool(processes=available_threads) as pool:
        results_list = pool.starmap(process_file, [(file_path, num_affinities) for file_path in file_paths])
    
    # Prepare to write results
    with open(f'{output_directory}/affinity_results.txt', 'w') as f:
        # Write header
        f.write("ligand protein affinities grid_size grid_space exhaustiveness random_seed\n")

        for file_name, affinities, grid_size, grid_space, exhaustiveness, random_seed in results_list:
            # Process file name to extract ligand and protein
            ligand = file_name.replace('.out', '').split('_')[0]
            protein = '_'.join(file_name.replace('.out', '').split('_')[1:])

            affinities_str = ",".join(map(str, affinities))  # Join with commas
            f.write(f"{ligand} {protein} {affinities_str} {grid_size} {grid_space} {exhaustiveness} {random_seed}\n")  # Add newline
            
    print(f"Affinity results have been saved to '{output_directory}/affinity_results.txt'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze affinity values from multiple files")
    parser.add_argument('--range', type=str, help="Range of folders in the format 'start-end'")
    parser.add_argument('--num', type=int, default=1, help="Number of affinities to extract (default is 1)")
    parser.add_argument('--threads', type=int, help="Number of threads to use (optional)")
    
    args = parser.parse_args()
    main(folder_range=args.range, num_affinities=args.num, num_threads=args.threads)
