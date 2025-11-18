import os

# 初始化一个空列表来保存满足条件的行
results = []

# 遍历指定路径下的所有txt文件
for filename in os.listdir('./out_8806'):
    if filename.endswith('.txt'):
        # 初始化两个列表分别用于保存num1和num2符合条件的行及其值
        num1_lines_values = []
        num2_lines_values = []
        
        with open(os.path.join('./out_8806', filename), 'r') as f:
            for line in f:
                # 根据'\t'进行文本切割
                parts = line.split('\t')
                # 确保有足够的分割元素
                if len(parts) >= 4:
                    try:
                        num1 = float(parts[2])  # 将字符串转换为浮点数
                        num2 = float(parts[3])
                    except ValueError:
                        continue  # 如果转换失败，则忽略这一行
                    
                    # 分别处理num1和num2，记录下符合条件的行和值
                    if num1 > 0.5:
                        num1_lines_values.append((line, num1))
                    if num2 > 0.5:
                        num2_lines_values.append((line, num2))
        
        # 分别对num1和num2的值进行排序，并取前5个 #这里实际上是取第一个
        num1_lines_values.sort(key=lambda x: x[1], reverse=True)
        num2_lines_values.sort(key=lambda x: x[1], reverse=True)
        top5_num1_lines = [line for line, value in num1_lines_values[:1]]
        top5_num2_lines = [line for line, value in num2_lines_values[:1] if line not in top5_num1_lines][:1]
        
        # 合并num1和num2的结果，避免重复
        file_results = top5_num1_lines + top5_num2_lines
        
        # 将当前文件的结果添加到总结果列表中
        results.extend(file_results)

# 把所有满足条件的行内容保存到新的本地文本中
with open('pocket_template_tmscore.txt', 'w') as f:
    for line in results:
        f.write(line)

