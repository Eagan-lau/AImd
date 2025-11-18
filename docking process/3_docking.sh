#!/bin/bash
#SBATCH -p big,low,smp01,smp02  ##节点名称
#SBATCH -J vina7   ##任务名称
#SBATCH -N 1  ##节点
#SBATCH -n 16  ##线程
#SBATCH -o vina-7.o  ##输出信息
#SBATCH -e vina-7.e  ##报错信息输出

python 4_reserve_docking.py -num 645 --exhaustiveness 16 --log_dir /public/agis/yanjianbin_group/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/4_docking_results/file_645
