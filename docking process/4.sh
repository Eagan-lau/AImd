#!/bin/bash
#SBATCH -p big,low,smp01,smp02  ##节点名称
#SBATCH -J results   ##vina
#SBATCH -N 1  ##节点
#SBATCH -n 16  ##线程
#SBATCH -o results.o  ##输出信息
#SBATCH -e results.e   ##报错信息输出


export PATH
 if [ -f "/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/etc/profile.d/conda.sh" ]; then
    . "/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/etc/profile.d/conda.sh"
 else
      export PATH="/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/bin:$PATH"
 fi
 conda activate docking   ##修改conda 环境名称
 export ALIGNMENT_TOOL_PATH=/public/agis/yanjianbin_group/liuyugeng/Anaconda/install/envs/docking/bin   ##修改conda 环境名称

python3 5_stastics_results.py --range 1-200 --num 9









