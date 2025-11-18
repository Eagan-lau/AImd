#!/bin/bash
#SBATCH -p low,big,smp01,smp02  ##节点名称
#SBATCH -J pdbqt1   ##任务名称
#SBATCH -N 1  ##节点
#SBATCH -n 10  ##线程
#SBATCH -o pdbqt-1.o  ##输出信息
#SBATCH -e pdbqt-1.e   ##报错信息输出


export PATH
 if [ -f "/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/etc/profile.d/conda.sh" ]; then
    . "/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/etc/profile.d/conda.sh"
 else
      export PATH="/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/bin:$PATH"
 fi
 conda activate docking   ##修改conda 环境名称
 export ALIGNMENT_TOOL_PATH=/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/envs/docking/bin   ##修改conda 环境名称

for num in {1..20}
do
   # 执行命令，使用$num替换-num参数
   python 1_PDBQT_process.py -mol protein -num $num
done











