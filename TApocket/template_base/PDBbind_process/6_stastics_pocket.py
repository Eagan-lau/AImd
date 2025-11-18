import os
import numpy as np

# 确保PyMOL的无头模式被启用
import pymol
from pymol import cmd  # 明确导入cmd

pymol.pymol_argv = ['pymol', '-qc']  # -q for quiet, -c for no GUI
pymol.finish_launching()

def extract_ca_coordinates(object_name):
    """提取指定对象的α-碳原子坐标"""
    model = cmd.get_model(f"{object_name} and name CA")
    ca_coords = np.array([[atom.coord[0], atom.coord[1], atom.coord[2]] for atom in model.atom])
    return ca_coords

def count_overlapping_atoms(coords1, coords2, threshold=2.0):
    """计算两组坐标中距离在阈值范围内的点的数量"""
    overlaps = 0
    for coord1 in coords1:
        distances = np.sqrt(np.sum((coords2 - coord1)**2, axis=1))
        overlaps += np.sum(distances <= threshold)
    return overlaps

pse_files_path = '/media/yugengliu/DATA1/pocket/8806/align_super'  # 指定.pse文件所在的路径
results_file_path = 'overlapping_percentage_results.txt'  # 结果文件路径

# 打开文件准备写入结果
with open(results_file_path, 'w') as results_file:
    for filename in os.listdir(pse_files_path):
        if filename.endswith('.pse'):
            file_path = os.path.join(pse_files_path, filename)
            
            pymol.cmd.load(file_path)  # 使用PyMOL打开.pse文件
            
            coords_pocket = extract_ca_coordinates('pocket')
            coords_one = extract_ca_coordinates('one')
            
            overlaps = count_overlapping_atoms(coords_pocket, coords_one)
            total_pocket_atoms = len(coords_pocket)
            overlap_percentage = (overlaps / total_pocket_atoms) * 100 if total_pocket_atoms > 0 else 0
            
            # 将结果写入文件
            results_file.write(f'{filename}: Overlapping α-carbon atoms percentage between "pocket" and "one": {overlap_percentage:.2f}%\n')
            
            pymol.cmd.delete('all')  # 清除当前PyMOL会话以加载下一个文件

# PyMOL会话结束
pymol.cmd.quit()

