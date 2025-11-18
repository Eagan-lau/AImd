#!/bin/bash
#SBATCH -p big,low,smp01,smp02  ## 节点队列
#SBATCH -J gate1_13   ## 任务名称
#SBATCH -N 1             ## 节点数
#SBATCH -n 2            ## 线程数
#SBATCH -o gate-1-13.o  ## 标准输出
#SBATCH -e gate-1-13.e  ## 错误输出

python 6_20250810_score.py --file-range 157
