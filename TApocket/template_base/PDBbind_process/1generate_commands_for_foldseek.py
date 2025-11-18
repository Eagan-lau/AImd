import os

# 设置源文件夹路径
source_dir = '/media/yugengliu/DATA1/pocket/8806_protein'
out_dir = '/home/yugengliu/DATA1/foldseek_databse/PDB_out_8806'
# 设置目标数据库路径
target_db = '/home/yugengliu/DATA1/foldseek_databse/PDB/pdb'


# 要写入的新文本文件名
commands_file = 'generated_commands_8806_PDB.txt'

# 创建一个空列表来保存所有命令
commands = []

# 遍历源文件夹中的所有文件
for filename in os.listdir(source_dir):
    if filename.endswith('.pdb'):  # 检查文件扩展名是否为.pdb
        # 从文件名中移除.pdb扩展名来获取{name}
        name = filename[:-4]
        
        # 格式化命令字符串
        command = f'foldseek easy-search {source_dir}/{filename} {target_db} {out_dir}/{name}.txt tmpFolder --format-output "query,target,qtmscore,ttmscore,alntmscore" --max-seqs 100'
        #command = f'foldseek easy-search {source_dir}/{filename} {target_db} {out_dir}/{name}.txt tmpFolder --format-output "query,target,qtmscore,ttmscore,alntmscore"'
        # 将命令添加到列表中
        commands.append(command)

# 将所有命令写入新文本文件
with open(commands_file, 'w') as file:
    for command in commands:
        file.write(command + '\n')  # 每个命令后面加上换行符

print(f'Commands have been written to {commands_file}')

