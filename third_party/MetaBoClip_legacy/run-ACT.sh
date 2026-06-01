#!/bin/bash
#SBATCH -p AMD_9A14  ##节点名称
#SBATCH -J distance_angle_act   ##任务名称
#SBATCH -N 1  ##节点
#SBATCH -n 16  ##线程
#SBATCH -o distance_angle_act.o  ##输出信息
#SBATCH -e distance_angle_act.e   ##报错信息输出

python run_family.py \
  --family act \
  --protein-dir /vol2/share/yanjianbin/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/1_PDBQT/protein/file_1 \
  --docking-dir /vol2/share/yanjianbin/liuyugeng/taxus_alphaflow_docking/ACT/yixianji/data_output/4_docking_results/file_1 \
  --out-dir /vol2/share/yanjianbin/liuyugeng/taxus_alphaflow_docking/docking-score-test-new/ACT-out