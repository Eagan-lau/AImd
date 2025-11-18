#!/bin/bash
#SBATCH -p big,low,smp01,smp02  ##节点名称
#SBATCH -J cavity1   ##任务名称
#SBATCH -N 1  ##节点
#SBATCH -n 16  ##线程
#SBATCH -o cavity1.o  ##输出信息
#SBATCH -e cavity1.e   ##报错信息输出


export PATH
 if [ -f "/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/etc/profile.d/conda.sh" ]; then
    . "/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/etc/profile.d/conda.sh"
 else
      export PATH="/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/bin:$PATH"
 fi
 conda activate docking   ##修改conda 环境名称
 export ALIGNMENT_TOOL_PATH=/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/envs/docking/bin   ##修改conda 环境名称

python 2_cavity_detected.py -num 3112 --threads 14 --timeout 2000









