import os
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from sklearn.cluster import KMeans

class ResidueSelect(Select):
    """自定义选择器，用于选择和保存残基到PDB文件"""
    def __init__(self, residues):
        self.residues = residues
    
    def accept_residue(self, residue):
        return residue in self.residues

def process_pdb_file(pdb_filename, output_directory):
    """处理单个PDB文件，并保存结果到指定目录"""
    pdb_parser = PDBParser(QUIET=True)
    structure = pdb_parser.get_structure('protein_structure', pdb_filename)

    # 检查是否有模型可用
    models = list(structure.get_models())
    if not models:
        print(f"No models found in {pdb_filename}. Skipping...")
        return
    
    model = models[0]  # 使用第一个模型

    alpha_carbon_coordinates = []
    residues = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                alpha_carbon = residue['CA']
                alpha_carbon_coordinates.append(alpha_carbon.get_coord())
                residues.append(residue)

    coords_array = np.array(alpha_carbon_coordinates)
    if len(coords_array) < 2:
        print(f"Skipping {pdb_filename} due to insufficient number of samples.")
        return

    # 使用K-means聚类，将数据分为两个类
    kmeans = KMeans(n_clusters=2)  #指定cluster的数量
    kmeans.fit(coords_array)
    cluster_labels = kmeans.labels_
    
    cluster_indices = set(cluster_labels)
    for cluster_label in cluster_indices:
        cluster_residues = [residues[i] for i, label in enumerate(cluster_labels) if label == cluster_label]
        
        pdbio = PDBIO()
        pdbio.set_structure(structure)
        cluster_filename = os.path.join(output_directory, f'{os.path.splitext(os.path.basename(pdb_filename))[0]}_cluster_{cluster_label+1}.pdb')
        pdbio.save(cluster_filename, ResidueSelect(cluster_residues))
        print(f'Saved cluster {cluster_label+1} to {cluster_filename}')

input_directory = './3_extracted_residues'  # 替换为实际的PDB文件目录路径
output_directory = './4_pocketminer_cluster'  # 替换为希望保存新PDB文件的目录路径

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

pdb_files = [f for f in os.listdir(input_directory) if f.endswith('.pdb')]
for pdb_file in pdb_files:
    process_pdb_file(os.path.join(input_directory, pdb_file), output_directory)
