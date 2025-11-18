import os
from Bio.PDB import PDBParser, PDBIO, Select

class BFactorSelect(Select):
    """自定义选择器，用于选择b-factor大于70的残基"""
    def accept_atom(self, atom):
        return atom.get_bfactor() > 70

def process_pdb_files(input_directory, output_directory):
    """处理输入路径下的所有PDB文件，并保存b-factor大于70的残基"""
    # 确保输出目录存在
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # 遍历输入目录下的所有PDB文件
    for pdb_file in os.listdir(input_directory):
        if pdb_file.endswith(".pdb"):
            input_file_path = os.path.join(input_directory, pdb_file)
            output_file_path = os.path.join(output_directory, pdb_file)

            # 解析PDB文件
            pdb_parser = PDBParser()
            structure = pdb_parser.get_structure('structure', input_file_path)

            # 保存b-factor大于70的残基到新的PDB文件
            pdbio = PDBIO()
            pdbio.set_structure(structure)
            pdbio.save(output_file_path, BFactorSelect())

            print(f"Processed {pdb_file}")

# 指定输入和输出路径
input_directory = './2_remove_disorder'
output_directory = './3_extracted_residues'

# 处理文件
process_pdb_files(input_directory, output_directory)
