#!/bin/bash
#SBATCH -p big,low,smp01,smp02  ## 节点队列
#SBATCH -J gate201_209   ## 任务名称
#SBATCH -N 1             ## 节点数
#SBATCH -n 2            ## 线程数
#SBATCH -o gate-201-209.o  ## 标准输出
#SBATCH -e gate-201-209.e  ## 错误输出

python 6_faxianliang.py --file-range 201-209
