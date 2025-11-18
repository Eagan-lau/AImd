#!/bin/bash
#SBATCH -p big,low,smp01,smp02  ##节点名称
#SBATCH -J 3   ##任务名称
#SBATCH -N 1  ##节点
#SBATCH -n 16  ##线程
#SBATCH -o 3.o  ##输出信息
#SBATCH -e 3.e  ##报错信息输出

python3 7_score_UGT-20250901.py
