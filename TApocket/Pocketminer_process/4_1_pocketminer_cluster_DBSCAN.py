import os
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

class ResidueSelect(Select):
    """自定义选择器，用于选择和保存残基到PDB文件"""
    def __init__(self, residues):
        self.residues = residues
    
    def accept_residue(self, residue):
        return residue in self.residues

def estimate_eps(coords_array, n_neighbors=10):
    """基于k-最近邻的平均距离来估计eps"""
    # 确保n_neighbors不大于样本数
    n_samples = coords_array.shape[0]
    if n_samples <= 1:
        # 对于只有一个样本的情况，返回一个默认的小eps值
        return 0.1
    n_neighbors = min(n_neighbors, max(2, n_samples // 2))  # 至少为2，且不超过样本数的一半

    nn = NearestNeighbors(n_neighbors=n_neighbors)
    nn.fit(coords_array)
    distances, _ = nn.kneighbors(coords_array)
    # 使用第n个最近邻距离的平均值作为eps，确保eps大于0
    eps = np.mean(distances[:, n_neighbors-1])
    return max(eps, 0.1)  # 确保eps至少为0.1，避免eps为0的情况

def process_pdb_file(pdb_filename, output_directory, eps_log_file):
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

    eps_value = estimate_eps(coords_array, n_neighbors=min(10, len(coords_array) // 2))
    min_samples = min(10, len(coords_array) // 2)
    pdb_basename = os.path.splitext(os.path.basename(pdb_filename))[0]

    with open(eps_log_file, 'a') as log:
        log.write(f"{pdb_basename}: eps={eps_value}, min_samples={min_samples}\n")

    dbscan = DBSCAN(eps=eps_value, min_samples=min(10, len(coords_array) // 2))
    dbscan.fit(coords_array)

    cluster_labels = set(dbscan.labels_)
    for cluster_label in cluster_labels:
        if cluster_label == -1:
            continue

        cluster_residues = [residues[i] for i, label in enumerate(dbscan.labels_) if label == cluster_label]
        
        pdbio = PDBIO()
        pdbio.set_structure(structure)
        cluster_filename = os.path.join(output_directory, f'{pdb_basename}_pocket_{cluster_label+1}.pdb')
        pdbio.save(cluster_filename, ResidueSelect(cluster_residues))
        print(f'Saved cluster {cluster_label+1} to {cluster_filename}')


pdb_directory = './3_extracted_residues'  # 替换为实际的PDB文件目录路径
output_directory = './4_pocketminer_cluster'  # 替换为希望保存新PDB文件的目录路径
output_directory2 = './'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

eps_log_file = os.path.join(output_directory2, 'eps_values.log')
open(eps_log_file, 'w').close()

pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]
for pdb_file in pdb_files:
    process_pdb_file(os.path.join(pdb_directory, pdb_file), output_directory, eps_log_file)
