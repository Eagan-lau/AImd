import os
from Bio.PDB import PDBParser, PDBIO, Select, DSSP

class DisorderSelect(Select):
    """这个类用于删除disorder区域的残基，如果这些残基连续且超过20个"""
    def __init__(self, disorder_regions):
        self.disorder_regions = disorder_regions

    def accept_residue(self, residue):
        # 如果残基在disorder区域内，则不接受（即删除）
        for start, end in self.disorder_regions:
            if start <= residue.id[1] <= end:
                return False
        return True

def process_pdb_files(input_dir, output_dir):
    parser = PDBParser()
    io = PDBIO()

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith(".pdb"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)

            # 读取并分析PDB文件
            structure = parser.get_structure(filename, input_path)
            model = structure[0]
            dssp = DSSP(model, input_path)

            # 查找disorder区域
            disorder_regions = []
            current_start = None
            current_count = 0
            for residue in model.get_residues():
                dssp_key = (residue.get_parent().id, residue.id)
                if dssp_key in dssp:
                    if dssp[dssp_key][2] == '-':
                        if current_start is None:
                            current_start = residue.id[1]
                        current_count += 1
                    else:
                        if current_count > 10:
                            disorder_regions.append((current_start, current_start + current_count - 1))
                        current_start = None
                        current_count = 0
                else:
                    if current_count > 10:
                        disorder_regions.append((current_start, current_start + current_count - 1))
                    current_start = None
                    current_count = 0

            # 如果最后的区域也是disorder区域且超过20个残基，则也需要添加
            if current_count > 10:
                disorder_regions.append((current_start, current_start + current_count - 1))

            # 删除disorder区域，并保存新的PDB文件
            selector = DisorderSelect(disorder_regions)
            io.set_structure(structure)
            io.save(output_path, selector)

if __name__ == "__main__":
    input_dir = "./1_b-factor/output_pdb"
    output_dir = "./2_remove_disorder"
    process_pdb_files(input_dir, output_dir)
