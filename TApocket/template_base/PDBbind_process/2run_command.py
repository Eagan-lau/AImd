import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

# 定义一个执行命令的函数，该函数将在线程池中的线程中被调用
def execute_command(command):
    try:
        # 使用subprocess.run执行命令，捕获stdout和stderr
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return f"Command executed successfully: {command}\n{result.stdout}"
    except subprocess.CalledProcessError as e:
        # 如果命令执行失败，返回错误信息
        return f"Error executing command: {command}\n{e.stderr}"

# 读取生成的命令文件
commands_file = 'generated_commands_8806_PDB.txt'
with open(commands_file, 'r') as file:
    commands = [command.strip() for command in file.readlines() if command.strip()]

# 设置最大线程数
max_threads = 28  # 例如，同时运行的最大线程数为10

# 使用ThreadPoolExecutor创建线程池
with ThreadPoolExecutor(max_workers=max_threads) as executor:
    # 将任务提交到线程池
    futures = [executor.submit(execute_command, command) for command in commands]
    
    # 等待所有任务完成，并打印结果
    for future in as_completed(futures):
        print(future.result())

print("All commands have been executed.")

