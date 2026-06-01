#!/bin/bash
#SBATCH -p AMD_9A14  ##节点名称
#SBATCH -J distance_angle_68   ##任务名称
#SBATCH -N 1  ##节点
#SBATCH -n 16  ##线程
#SBATCH -o distance_angle.o  ##输出信息
#SBATCH -e distance_angle.e   ##报错信息输出

python run_family.py \
  --family ugt \
  --protein-dir /vol2/share/yanjianbin/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/1_PDBQT/protein/file_1 \
  --docking-dir /vol2/share/yanjianbin/liuyugeng/taxus_alphaflow_docking/UDPGT/data_output/4_docking_results/file_1 \
  --out-dir /vol2/share/yanjianbin/liuyugeng/taxus_alphaflow_docking/docking-score-test-new/UDPGT-out