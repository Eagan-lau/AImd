#!/bin/bash
#SBATCH -p big,low,smp01,smp02  ## 节点队列
#SBATCH -J gate1_5   ## 任务名称
#SBATCH -N 1             ## 节点数
#SBATCH -n 2            ## 线程数
#SBATCH -o gate-1-5.o  ## 标准输出
#SBATCH -e gate-1-5.e  ## 错误输出

python 6_20250815_CH_only.py --file-range 631
