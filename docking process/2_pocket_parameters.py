"""
Author: ygliu
Description: This script is designed for generating coordinate information for protein-ligand docking. Specifically, it extracts the center coordinates of the pocket from each cavity in the protein structure. Then, it calculates the dimensions for docking and writes these parameters to a text file. It uses command-line arguments to specify details like the maximum number (30) of cavities to consider, the ligand name, and other parameters.
"""

import os
import argparse
from Bio.PDB import PDBParser
import math

# Function to get the list of input file names from the given directory
def get_input_file_names(path):
    names = []
    for filename in os.listdir(path):
        #if filename.endswith(".input"):
            #names.append(os.path.splitext(filename)[0])
        if filename.endswith(".pdb"):
            names.append(os.path.splitext(filename)[0].split('_cavity_')[0])
        elif filename.startswith('pocketminer_'):
            names.append(os.path.splitext(filename)[0].replace('pocketminer_',''))
    return names

# Function to define the XYZ coordinates of a pocket
def find_XYZ_coordinates(structure, center_x, center_y, center_z):
    # Initialize min and max coordinates as the center to ensure that the difference is calculated from the center
    min_x, max_x = center_x, center_x
    min_y, max_y = center_y, center_y
    min_z, max_z = center_z, center_z

    # Traverse through atoms to update min and max coordinates
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x, y, z = atom.coord
                    min_x = min(min_x, x)
                    max_x = max(max_x, x)
                    min_y = min(min_y, y)
                    max_y = max(max_y, y)
                    min_z = min(min_z, z)
                    max_z = max(max_z, z)

    # Calculate the maximum distance from the center for each axis
    max_dist_x = max(center_x - min_x, max_x - center_x)
    max_dist_y = max(center_y - min_y, max_y - center_y)
    max_dist_z = max(center_z - min_z, max_z - center_z)

    return max_dist_x, max_dist_y, max_dist_z

# Function to calculate the center of the pocket in the PDB structure
def calculate_pocket_center(structure):
    x = 0
    y = 0
    z = 0
    count = 0
    # Traverse through atoms and calculate pocket center
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    count += 1
                    x += atom.coord[0]
                    y += atom.coord[1]
                    z += atom.coord[2]
    if count == 0:  # To avoid division by zero
        return 0, 0, 0
    return x / count, y / count, z / count

# Main script starts here
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate coordinate information', add_help=False)
    
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    
    required.add_argument('-proteinfolder_num', type=int, required=True, help='Number for protein folder')
    required.add_argument('-ligand', type=str, required=True, help='Name of the ligand')
    optional.add_argument('--max_cavity', type=int, default=3, help='Max cavity index to be processed,default=3')
    optional.add_argument('--output_path_num', type=str, help='The output path for pocket parameters and docking results,default=proteinfolder_num')
    optional.add_argument('--ligandfolder_num', type=int, default=1, help='Uses the number to construct the ligand path, ./data_output/1_PDBQT/ligand/file_{num},default=1')
    optional.add_argument('--size_xyz', type=str, default=None, help='Size dimensions separated by commas,e.g. 15,15,15')
    optional.add_argument('--pocket_real_size', action='store_true', help='Use the actual size of the Cavity detection pocket. If not enabled, max dimension is 30 Angstroms. Advice: Increase "exhaustiveness" in docking process if enabled')
    optional.add_argument('--min_size', type=int, default=5, help='Minimum size for x, y, z dimensions, default=5')
    optional.add_argument('--margin', type=int, default=2, help='Additional margin for the pocket size, default=2')
    optional.add_argument("-h", "--help", action="help", help="Display this help message")
    args = parser.parse_args()
    
    if args.output_path_num is None:
        args.output_path_num = args.proteinfolder_num
        
    current_path = os.getcwd()
    pdbqt_protein_path = f'{current_path}/data_output/1_PDBQT/protein/file_{args.proteinfolder_num}'
    pdbqt_ligand_path = f'{current_path}/data_output/1_PDBQT/ligand/file_{args.ligandfolder_num}'
    cavity_path = f'{current_path}/data_output/2_cavity_output/cavity_{args.proteinfolder_num}'
    output_path = f'{current_path}/data_output/3_pocket_parameters/file_{args.output_path_num}'
    docking_output_path = f'{current_path}/data_output/4_docking_results/file_{args.output_path_num}'
    
    if not os.path.exists(docking_output_path):
        os.makedirs(docking_output_path)
        
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    input_names = get_input_file_names(cavity_path)
    # Initialize PDBParser from BioPython
    p = PDBParser()
    
    # Initialize a counter for structure ID
    structure_id_counter = 1

    # Iterate through each input file and calculate coordinate parameters
    for name in input_names:
        for i in range(1, args.max_cavity + 1):
            if os.path.exists(f"{cavity_path}/pocketminer_{name}.pdb"):
                cavity_file = f"{cavity_path}/pocketminer_{name}.pdb"
                new_name = f"{args.ligand}_pocketminer_{name}"
                i = '1'
            else:
                cavity_file = f"{cavity_path}/{name}_cavity_{i}.pdb"
                new_name = f"{args.ligand}_{name}_cavity_{i}"

            if not os.path.exists(cavity_file):
                continue

            # Use a simple and universal structure ID
            structure_id = f"structure_{structure_id_counter}"
            structure_id_counter += 1

            structure = p.get_structure(structure_id, cavity_file)
            
            center_x, center_y, center_z = calculate_pocket_center(structure)

            if args.size_xyz:
                size_x, size_y, size_z = map(int, args.size_xyz.split(','))
            else:
                if args.pocket_real_size:
                    max_dist_x, max_dist_y, max_dist_z = find_XYZ_coordinates(structure, center_x, center_y, center_z)
                    size_x = 2 * math.ceil(max_dist_x) + args.margin
                    size_y = 2 * math.ceil(max_dist_y) + args.margin
                    size_z = 2 * math.ceil(max_dist_z) + args.margin
                    
                else:
                    max_dist_x, max_dist_y, max_dist_z = find_XYZ_coordinates(structure, center_x, center_y, center_z)
                    size_x = min(2 * math.ceil(max_dist_x) + args.margin, 30)
                    size_y = min(2 * math.ceil(max_dist_y) + args.margin, 30)
                    size_z = min(2 * math.ceil(max_dist_z) + args.margin, 30)
            
                #Check and adjust to minimum size
                size_x = max(size_x, args.min_size)
                size_y = max(size_y, args.min_size)
                size_z = max(size_z, args.min_size)
                     
            lines = [
                f"receptor = {pdbqt_protein_path}/{name}.pdbqt\n",
                f"ligand = {pdbqt_ligand_path}/{args.ligand}.pdbqt\n",
                f"out = {docking_output_path}/{args.ligand}@{name}_c{i}_out.pdbqt\n",
                f"center_x = {round(center_x, 2)}\n",
                f"center_y = {round(center_y, 2)}\n",
                f"center_z = {round(center_z, 2)}\n",
                f"size_x = {size_x}\n",
                f"size_y = {size_y}\n",
                f"size_z = {size_z}\n"
            ]

            with open(f'{output_path}/{new_name}.txt', 'w') as f:
                f.writelines(lines)
        #print(f"Completed coordinate generation for: {name}")
    print("Coordinate information generation for all cavities completed.")