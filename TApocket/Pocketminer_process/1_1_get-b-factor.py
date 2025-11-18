from Bio.PDB import PDBParser, PDBIO
import os

def update_b_factors(pdb_path, b_factors, output_path):
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_path)

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # 确保我们不会越界
                    if atom.get_serial_number() - 1 < len(b_factors):
                        atom.set_bfactor(b_factors[atom.get_serial_number() - 1])

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path)

# 定义你的目录路径
input_txt_dir = './0_data_input/pocketminer_out_txt'
input_pdb_dir = './0_data_input/protein_pdb'
output_dir = './1_b-factor/output_pdb'

# 确保输出目录存在
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 遍历.txt文件
for txt_file in os.listdir(input_txt_dir):
    if txt_file.endswith('-predictions.txt'):
        base_name = txt_file.replace('-predictions.txt', '')
        pdb_name = base_name + '.pdb'
        b_factors_file = os.path.join(input_txt_dir, txt_file)
        pdb_file_path = os.path.join(input_pdb_dir, pdb_name)
        output_pdb_path = os.path.join(output_dir, pdb_name)
        
        # 读取b-factor值
        with open(b_factors_file, 'r') as file:
            b_factors = [float(line.strip()) * 100 for line in file.readlines()]
        
        # 更新PDB文件中的b-factor值
        update_b_factors(pdb_file_path, b_factors, output_pdb_path)

print("所有修改后的PDB文件已保存到" + output_dir)
