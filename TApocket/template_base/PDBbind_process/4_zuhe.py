# 导入必要的库
import os

# 初始化一个字典来保存protein和name的对应关系
protein_dict = {}

# 读取stastics_output.txt的内容，逐行进行分析
with open('pocket_template_tmscore.txt', 'r') as f:
    for line in f:
        # 先根据'\t'进行切割，切割后取第一个元素{protein}跟第二个元素{pocket}
        protein, pocket = line.split('\t')[0:2]
        # 第二个元素{pocket}再根据'.cif'进行切割后取第一个元素，然后这个元素加上'.pdb'构成完整的文件名{name}
        name = pocket.split('.pdb')[0] + '.pdb'
        # 生成1中{protein}跟{name}的对应关系，然后相同的{protein}合并在一起
        if protein not in protein_dict:
            protein_dict[protein] = [name]
        else:
            protein_dict[protein].append(name)

# 将这些处理后的新的内容保存到新的本地文本中
with open('1output.txt', 'w') as f:
    for protein, names in protein_dict.items():
        # 一个{protein}占据一行的空间，存储格式为{protein} {name} {name} {name}
        f.write(protein + ' ' + ' '.join(names) + '\n')

