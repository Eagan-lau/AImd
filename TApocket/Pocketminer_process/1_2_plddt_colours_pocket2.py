import os
import pymol
from pymol import cmd, stored

def process_pdb(pdb_file, output_folder):
    # 确保PyMOL以无头模式启动（无GUI）
    pymol.finish_launching(['pymol', '-qc'])  # '-q' for quiet, '-c' for command line
    
    # 构建输出文件路径
    base_name = os.path.basename(pdb_file)
    output_file = os.path.join(output_folder, base_name.replace('.pdb', '.pse'))

    # 加载PDB文件
    cmd.load(pdb_file, "structure")

    # 应用颜色设置
    cmd.spectrum("b", "blue_white_red", minimum=0, maximum=100, selection='structure')

    # 保存为PSE格式，包含原始结构
    cmd.save(output_file, "structure")

    # 清理
    cmd.delete("all")

if __name__ == "__main__":
    input_folder = './1_b-factor/output_pdb'  # 输入目录路径
    output_folder = './1_b-factor/output_color'  # 输出目录路径

    # 确保输出目录存在
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 遍历输入目录下的所有.pdb文件
    for file in os.listdir(input_folder):
        if file.endswith('.pdb'):
            pdb_file = os.path.join(input_folder, file)
            process_pdb(pdb_file, output_folder)
    print("处理完成。")
