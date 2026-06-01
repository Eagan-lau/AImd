# AImd 模块接口与 manifest 规范

本文件记录 AImd 当前模块之间的标准接口。后续替换工具时，应优先保持这些接口不变。

---

## 1. 总体数据流

```text
RGPC
  outputs:
    data/protein/protein_manifest.csv
    data/cluster/RGPC/clusters.tsv
    data/cluster/RGPC/representatives.tsv

TApocketBridge
  inputs:
    data/protein/protein_manifest.csv
  outputs:
    data/pocket/pocket_manifest.csv
    data/pocket/tapocket_run_manifest.csv

DockingHub broad
  inputs:
    data/protein/protein_manifest.csv
    data/pocket/pocket_manifest.csv
    data/ligand/ligand_manifest.csv
  outputs:
    data/docking_out/docking_result_manifest.csv

ClusterScore
  inputs:
    data/docking_out/docking_result_manifest.csv
  outputs:
    data/scoring/ClusterScore/clusterscore_results.xlsx
    data/scoring/ClusterScore/top10_clusters.csv

RefinementHub
  inputs:
    data/scoring/ClusterScore/top10_clusters.csv
    data/protein/protein_manifest.csv
  outputs:
    data/refinement/selected_protein_manifest.csv
    data/refinement/refined_docking.generated.yaml
    data/refined/docking_out/docking_result_manifest.csv

MetaBoClipBridge
  inputs:
    data/refined/docking_out/docking_result_manifest.csv
  outputs:
    data/metaboclip/results/metaboclip_final_ranking.csv
```

---

## 2. protein_manifest.csv

来源：RGPC  
下游：TApocketBridge、DockingHub、RefinementHub

建议字段：

```csv
protein_id,cluster_id,batch_id,source_pdb,protein_path,is_representative,file_action,status
```

字段说明：

| 字段 | 必需 | 说明 |
|---|---:|---|
| protein_id | 是 | 蛋白唯一 ID |
| cluster_id | 是 | RGPC 聚类 ID |
| batch_id | 是 | `file_1`, `file_2` 等 |
| source_pdb | 否 | 原始结构路径 |
| protein_path | 是 | 当前用于下游模块的结构路径 |
| is_representative | 否 | 是否 cluster representative |
| file_action | 否 | copy/symlink |
| status | 是 | success/failed 等 |

---

## 3. pocket_manifest.csv

来源：TApocketBridge  
下游：DockingHub

关键字段：

```csv
protein_id,cluster_id,batch_id,pocket_id,pocket_rank,center_x,center_y,center_z,size_x,size_y,size_z,final_score,protein_path,pocket_pdb_path,pocket_json_path,box_yaml_path,status
```

字段说明：

| 字段 | 必需 | 说明 |
|---|---:|---|
| protein_id | 是 | 必须与 protein_manifest.csv 对应 |
| cluster_id | 否 | 用于结果追踪 |
| batch_id | 是 | 必须与 protein 的 batch 对应 |
| pocket_id | 是 | pocket 唯一 ID |
| pocket_rank | 是 | pocket 排名，DockingHub 默认取 top 1 |
| center_x/y/z | 是 | Vina box 中心 |
| size_x/y/z | 是 | Vina box 尺寸 |
| final_score | 否 | 口袋预测得分 |
| pocket_pdb_path | 否 | pocket 结构文件 |
| box_yaml_path | 否 | box 配置 |
| status | 是 | success/failed |

---

## 4. ligand_manifest.csv

来源：分子处理模块或手动准备  
下游：DockingHub

建议字段：

```csv
ligand_id,batch_id,ligand_path,smiles,name,status
```

字段说明：

| 字段 | 必需 | 说明 |
|---|---:|---|
| ligand_id | 是 | ligand 唯一 ID |
| batch_id | 是 | ligand batch，一般可以全放 file_1 |
| ligand_path | 是 | PDBQT 路径 |
| smiles | 否 | 原始 SMILES |
| name | 否 | 分子名称 |
| status | 是 | success/failed |

---

## 5. docking_result_manifest.csv

来源：DockingHub  
下游：ClusterScore、MetaBoClipBridge

关键字段：

```csv
job_id,ligand_id,protein_id,cluster_id,batch_id,conformer_id,pocket_id,pocket_rank,receptor_pdbqt_path,ligand_pdbqt_path,config_path,out_pose_path,log_path,center_x,center_y,center_z,size_x,size_y,size_z,status,return_code,message,best_affinity,affinities,n_affinities,grid_size,grid_space,exhaustiveness,random_seed,pose_exists
```

字段说明：

| 字段 | 必需 | 说明 |
|---|---:|---|
| job_id | 是 | ligand/protein/pocket/conformer 组合 ID |
| ligand_id | 是 | ligand ID |
| protein_id | 是 | protein ID |
| cluster_id | 是 | cluster ID |
| batch_id | 是 | file_* |
| conformer_id | 是 | conf_0, conf_1 等 |
| pocket_id | 是 | pocket ID |
| receptor_pdbqt_path | 是 | receptor PDBQT |
| ligand_pdbqt_path | 是 | ligand PDBQT |
| out_pose_path | 是 | Vina pose 输出 |
| log_path | 是 | Vina `.out` 日志 |
| best_affinity | 是 | 第一个构象结合能 |
| affinities | 否 | 多个构象 affinity |
| center_x/y/z | 是 | box center |
| size_x/y/z | 是 | box size |
| status | 是 | success/skipped/failed/not_run |

---

## 6. ClusterScore 输出

`best_affinity_long.csv`：

```csv
job_id,ligand_id,protein_id,cluster_id,batch_id,pocket_id,conformer_id,best_affinity,...
```

`protein_ligand_matrix.csv`：

```csv
name,cluster,LIG0001,LIG0002,...
protein_A,C000001,-8.2,-5.1,...
```

`cluster_binding_statistics.csv`：

```csv
rank,cluster,alpha,beta,n_proteins,avg_strong_ratio,avg_weak_ratio,avg_none_ratio,total_strong_count,total_weak_count,total_none_count,total_binding_count,best_affinity_min,mean_affinity,composite_score,formula
```

`top10_clusters.csv`：

```csv
rank,cluster,composite_score,...
```

---

## 7. RefinementHub 输出

`selected_clusters.csv`：

```csv
selected_rank,cluster_id,rank,cluster,composite_score,...
```

`selected_protein_manifest.csv`：

```csv
protein_id,cluster_id,batch_id,protein_path,is_representative,status,selected_for_refinement,cluster_rank
```

`refined_docking.generated.yaml`：

由 `configs/Docking/docking.yaml` 与 `configs/Refinement/refine_from_clusterscore.yaml` 中的 `docking_overrides` 深度合并生成。

---

## 8. MetaBoClipBridge 输出

`metaboclip_staging_manifest.csv`：

```csv
family,ligand_id,protein_id,cluster_id,protein_name,pocket_id,conformer_id,job_id,best_affinity,protein_dir,docking_dir,staged_receptor_path,staged_pose_path,staged_log_path
```

最终结果：

```csv
metaboclip_final_ranking.csv
```

字段取决于 MetaBoClip 家族规则输出，通常包括：

```text
family
ligand_id
protein_name/protein_id
pose-level features
gating result
protein_score_norm
max_s_r
rank
```

---

## 9. 工具替换原则

替换某个模块内部工具时，应保证：

1. 输入 manifest 可读取；
2. 输出 manifest 字段不变；
3. 路径字段使用相对 AImd root 或绝对路径；
4. `status` 字段明确表示成功或失败；
5. 不成功的对象不要进入下游计算。

---
