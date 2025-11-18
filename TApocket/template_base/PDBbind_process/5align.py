import os
import pymol2

# 文件和目录路径
input_file = '1output.txt'
protein_dir = './8714_8806_protein_PDBbind'
output_dir = './align_super'

# 确保输出目录存在
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 使用 PyMOL 进行蛋白质对齐
def align_proteins(one, protein, pocket, output_file):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(os.path.join(protein_dir, f"{one}"), 'one')
        pymol.cmd.load(os.path.join(protein_dir, f"{protein}"), 'protein')
        pymol.cmd.load(os.path.join(protein_dir, f"{pocket}.pdb"), 'pocket')

        # 对齐操作：以 'one' 为 mobile，以 'protein' 为 target，'pocket' 保持不动
        pymol.cmd.super('one', 'protein')

        # 保存包含所有三个蛋白质的会话
        pymol.cmd.save(output_file)

# 读取并处理文件
with open(input_file, 'r') as f:
    for line in f:
        elements = line.strip().split(' ')
        one = elements[0]
        two_base = elements[1].split('.pdb')[0]  # 移除'.pdb'后缀得到基础名称
        protein = two_base + '.pdb'  # 第二个蛋白的完整名称

        # TODO: 需要一个明确的规则来定义如何从'two_base'获取'pocket'名称
        # 这里假设'pocket'名称与'two_base'相同
        pocket = two_base.replace('protein','pocket')

        # 构建输出文件名
        output_file = os.path.join(output_dir, f"{one}_{protein}.pse")

        # 对齐蛋白质并保存结果
        align_proteins(one, protein, pocket, output_file)

print("Alignment completed and saved.")

