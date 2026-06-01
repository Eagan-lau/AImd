# AImd 使用手册

版本：AImd integrated workflow  
核心流程：RGPC → TApocketBridge → DockingHub → ClusterScore → RefinementHub → refined DockingHub → MetaBoClipBridge

---

## 0. 总体说明

AImd 是一个以分子对接为核心枢纽的模块化工程包。当前版本已经把以下模块串联起来：

```text
蛋白结构库
  ↓
RGPC：Foldseek + HipMCL 结构聚类
  ↓
protein/file_* 和 protein_manifest.csv
  ↓
TApocketBridge：基于原始蛋白结构预测 pocket
  ↓
pocket_manifest.csv
  ↓
DockingHub：粗略 docking
  ↓
docking_result_manifest.csv
  ↓
ClusterScore：提取最优结合能并排序重要 cluster
  ↓
RefinementHub：筛选 Top cluster
  ↓
DockingHub：精细 docking，支持 AlphaFlow ensemble 和 cofactor 映射
  ↓
MetaBoClipBridge：构象门控、家族特征打分、pose/protein 排序
```

模块之间通过标准 manifest 文件连接。后续修改某个模块时，只要保持对应 manifest 的字段不变，下游模块通常不需要改。

---

## 1. 工程目录结构

解压后进入：

```bash
cd AImd
```

主要目录如下：

```text
AImd/
├── RGPC/                       # 蛋白结构聚类模块
├── TApocket/                   # TApocket 原始工程包
├── TApocketBridge/             # AImd 与 TApocket 的批处理接口
├── DockingHub/                 # ensemble/cofactor-aware docking 模块
├── ScoringHub/                 # ClusterScore 模块
├── RefinementHub/              # 根据 ClusterScore 筛选 top cluster 并复对接
├── MetaBoClip/                 # MetaBoClip 原始工程包
├── MetaBoClipBridge/           # AImd 与 MetaBoClip 的接口
├── orchestrator/               # 串联式工作流
├── registry/                   # 预留工具热插拔接口说明
├── configs/                    # 所有模块配置文件
├── data/                       # 输入、中间结果、最终结果
├── docs/                       # 使用手册和接口说明
├── run_rgpc.py
├── run_tapocket_batch.py
├── run_docking.py
├── run_clusterscore.py
├── run_refinement.py
├── run_metaboclip_bridge.py
└── run_full_iterative_metaboclip.py
```

---

## 2. 环境准备

### 2.1 Python 环境

建议使用 conda 或 venv。

```bash
cd AImd

conda create -n aimd python=3.10 -y
conda activate aimd

pip install -r requirements.txt
```

当前 `requirements.txt` 主要包含：

```text
PyYAML
biopython
pandas
numpy
matplotlib
```

如果你要运行 MetaBoClip 的真实几何分析，还需要保证 PyMOL 或 pymol2 可用。不同服务器环境下 PyMOL 的安装方式不同，建议先确认：

```bash
pymol -cq
```

如果命令不可用，需要先安装或把 PyMOL 加入 PATH。

### 2.2 外部工具

完整流程可能需要以下外部工具：

| 工具 | 主要使用模块 | 用途 |
|---|---|---|
| `foldseek` | RGPC、DockingHub cofactor | 结构检索、结构相似性计算、cofactor 模板选择 |
| `hipmcl` | RGPC | 基于结构相似图的聚类 |
| `pymol` | TApocket、DockingHub、MetaBoClip | 结构对齐、cofactor 映射、构象几何分析 |
| `prepare_receptor` | DockingHub | receptor PDB/PDBQT 准备 |
| `vina` | DockingHub | AutoDock Vina 对接 |
| `alphaflow` | DockingHub ensemble，可选 | 蛋白多构象采样 |

可以用下面命令检查工程布局和常用外部工具：

```bash
python validate_aimd_layout.py --root .
```

如果希望缺少外部工具时直接报错：

```bash
python validate_aimd_layout.py --root . --strict-tools
```

---

## 3. 必须遵守的运行规则

### 3.1 运行命令前先进入 AImd 根目录

所有默认配置中的：

```yaml
aimd_root: .
```

都假设你当前工作目录是 `AImd/`。

因此请始终：

```bash
cd AImd
```

再运行命令。

### 3.2 口袋预测与 ensemble 的坐标规则

AImd 当前采用固定规则：

```text
无论是否启用 ensemble docking，都先用原始蛋白结构预测 pocket。
AlphaFlow 或其他工具产生的 ensemble conformer 不重新预测 pocket。
每个 ensemble conformer 必须对齐回原始蛋白结构坐标系。
对接时复用原始 TApocket pocket box。
```

这样可以保证：

```text
pocket center / box size 始终位于原始结构坐标系；
ensemble conformer 经过对齐后与原始 pocket 坐标一致；
不会因为 AlphaFlow 输出坐标系漂移导致 docking box 偏移。
```

默认配置中已经设置：

```yaml
ensemble:
  align_to_reference:
    enabled: true
    fallback_to_unaligned: false
```

生产运行建议保持 `fallback_to_unaligned: false`。如果 PyMOL 对齐失败，则该 conformer 不应进入 docking。

---

## 4. 数据准备总览

完整流程需要以下输入数据：

```text
AImd/data/protein_structure/cleaned_pdb/     # 原始蛋白结构库
AImd/data/ligand/ligand_manifest.csv         # ligand manifest
AImd/data/ligand/file_*/                     # ligand PDBQT
AImd/data/cofactor/file_*/                   # 可选，cofactor 模板
AImd/TApocket/database_process/...           # 可选，TApocket 模板库和 M-CSA 数据库
```

其中蛋白结构库是当前流程的起点。分子处理模块后续可以单独接入，但在运行 docking 前必须已经有 ligand PDBQT 和 ligand manifest。

---

## 5. Step 1：放入蛋白结构库

### 5.1 输入位置

将蛋白结构放入：

```text
AImd/data/protein_structure/cleaned_pdb/
```

示例：

```text
AImd/data/protein_structure/cleaned_pdb/
├── TcChr01b00663.t1.pdb
├── TcChr01b00777.t1.pdb
├── TcChr05a01234.t1.pdb
└── ...
```

支持扩展名由 `configs/RGPC/rgpc.yaml` 控制：

```yaml
input:
  file_extensions: [".pdb", ".cif", ".mmcif"]
```

### 5.2 可选：使用 structure_manifest.csv

如果你想显式指定 protein_id 和结构路径，可以创建：

```text
AImd/data/protein_structure/structure_manifest.csv
```

格式：

```csv
protein_id,pdb_path,status
TcChr01b00663.t1,data/protein_structure/cleaned_pdb/TcChr01b00663.t1.pdb,success
TcChr01b00777.t1,data/protein_structure/cleaned_pdb/TcChr01b00777.t1.pdb,success
```

如果这个文件存在，RGPC 会优先读取它；否则扫描 `cleaned_pdb/`。

---

## 6. Step 2：运行 RGPC 结构聚类

### 6.1 配置文件

```text
AImd/configs/RGPC/rgpc.yaml
```

关键参数：

```yaml
foldseek:
  executable: foldseek
  mode: all_vs_all
  threads: 32

filter:
  min_qtmscore: 0.5
  min_ttmscore: 0.5
  weight_method: min_qtmscore_ttmscore

hipmcl:
  executable: hipmcl
  inflation: 2.0

batching:
  output_to_protein_dir: true
  max_proteins_per_batch: 500
  keep_cluster_together: true
  file_action: copy
```

### 6.2 运行命令

```bash
python run_rgpc.py --config configs/RGPC/rgpc.yaml
```

### 6.3 主要输出

RGPC 聚类结果：

```text
AImd/data/cluster/RGPC/
├── rgpc_input_manifest.csv
├── foldseek_search/all_vs_all.tsv
├── graph/structure_edges.tsv
├── graph/structure_edges.abc
├── hipmcl/hipmcl_output.txt
├── clusters.tsv
├── representatives.tsv
└── cluster_summary.csv
```

下游使用的 protein 文件：

```text
AImd/data/protein/
├── file_1/
│   ├── protein_A.pdb
│   └── protein_B.pdb
├── file_2/
│   └── ...
└── protein_manifest.csv
```

`protein_manifest.csv` 是后续 TApocket 和 DockingHub 的关键接口，格式：

```csv
protein_id,cluster_id,batch_id,source_pdb,protein_path,is_representative,file_action,status
TcChr01b00663.t1,C000001,file_1,/abs/path/TcChr01b00663.t1.pdb,data/protein/file_1/TcChr01b00663.t1.pdb,true,copy,success
```

### 6.4 后续替换其它聚类工具

RGPC 当前使用：

```text
Foldseek + HipMCL
```

后续可替换为：

```text
Foldseek easy-cluster
MMseqs2 + MCL
TM-align/US-align + Leiden
自定义结构相似性网络
```

只需要保证最后输出：

```text
data/protein/protein_manifest.csv
data/cluster/RGPC/clusters.tsv
data/cluster/RGPC/representatives.tsv
```

即可继续接入后续模块。

---

## 7. Step 3：运行 TApocket 口袋预测

### 7.1 口袋预测输入

TApocketBridge 读取：

```text
AImd/data/protein/protein_manifest.csv
```

并逐个处理其中的 `protein_path`。

### 7.2 TApocket 模板库和数据库放置位置

如果你使用 TApocket 的模板检索、M-CSA 过滤或 AI pocket 预测，需要根据 TApocket 内部配置准备数据。当前默认位置包括：

```text
AImd/TApocket/database_process/PLDB_redundancy/template_merged/
AImd/TApocket/database_process/MCSA/reference_active_sites/
AImd/TApocket/models/deeppocket_db/
AImd/TApocket/indexes/
```

这些目录中默认只有占位文件。真实运行前需要放入你的模板结构、M-CSA 活性位点或模型文件，并根据：

```text
AImd/TApocket/configs/tapocket_template_v1.yaml
```

调整具体路径。

### 7.3 配置文件

```text
AImd/configs/TApocket/tapocket_batch.yaml
```

关键参数：

```yaml
paths:
  protein_manifest: data/protein/protein_manifest.csv
  tapocket_project_dir: TApocket
  tapocket_config: TApocket/configs/tapocket_template_v1.yaml
  out_dir: data/pocket

selection:
  mode: all

execution:
  workers: 1
  continue_on_error: true
```

如果只想对代表蛋白预测 pocket：

```yaml
selection:
  mode: representatives_only
```

### 7.4 运行命令

```bash
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml
```

### 7.5 主要输出

```text
AImd/data/pocket/
├── file_1/
│   ├── protein_A/
│   │   ├── final_pockets.pdb
│   │   ├── final_pockets.json
│   │   ├── final_pocket_residues.tsv
│   │   ├── final_boxes.tsv
│   │   ├── summary.tsv
│   │   └── run_summary.json
│   └── ...
├── tapocket_run_manifest.csv
├── pocket_manifest.csv
├── pocket_batch_summary.csv
└── tapocket_batch_report.json
```

`pocket_manifest.csv` 是 docking 的关键输入：

```csv
protein_id,cluster_id,batch_id,pocket_id,pocket_rank,center_x,center_y,center_z,size_x,size_y,size_z,final_score,protein_path,pocket_pdb_path,pocket_json_path,box_yaml_path,status
```

### 7.6 后续替换其它口袋预测工具

可以替换为：

```text
P2Rank
fpocket
PocketMiner
DeepPocket
自研 TApocket 新版本
```

只要新模块最后写出：

```text
data/pocket/pocket_manifest.csv
```

并包含：

```text
protein_id, batch_id, pocket_id, pocket_rank,
center_x, center_y, center_z,
size_x, size_y, size_z, status
```

DockingHub 就可以继续运行。

---

## 8. Step 4：准备 ligand 数据

当前包尚未把分子处理模块接入主流程，因此 docking 前需要你手动或通过已有脚本准备 ligand PDBQT。

### 8.1 ligand 文件夹

建议按 batch 放入：

```text
AImd/data/ligand/
├── file_1/
│   ├── LIG0001.pdbqt
│   ├── LIG0002.pdbqt
│   └── ...
├── file_2/
│   └── ...
└── ligand_manifest.csv
```

### 8.2 ligand_manifest.csv 格式

```csv
ligand_id,batch_id,ligand_path,smiles,name,status
LIG0001,file_1,data/ligand/file_1/LIG0001.pdbqt,CCO,ligand_1,success
LIG0002,file_1,data/ligand/file_1/LIG0002.pdbqt,CCC,ligand_2,success
```

至少需要：

```text
ligand_id
batch_id
ligand_path
status
```

DockingHub 也能在没有 `ligand_manifest.csv` 时扫描：

```text
data/ligand/file_*/*.pdbqt
```

但为了保证可追踪性，建议始终提供 manifest。

---

## 9. Step 5：粗略分子对接 DockingHub

### 9.1 输入文件

粗略 docking 需要：

```text
data/protein/protein_manifest.csv
data/pocket/pocket_manifest.csv
data/ligand/ligand_manifest.csv
```

### 9.2 配置文件

```text
AImd/configs/Docking/docking.yaml
```

粗略 docking 默认：

```yaml
ensemble:
  enabled: false

cofactor:
  enabled: false

docking:
  engine: vina
  run: true
```

如果你只想生成 docking task 和 config，不真实运行 Vina：

```yaml
docking:
  run: false
```

### 9.3 运行命令

```bash
python run_docking.py --config configs/Docking/docking.yaml
```

### 9.4 主要输出

```text
AImd/data/docking_configs/file_*/*.txt
AImd/data/docking_tasks/docking_task_manifest.csv
AImd/data/docking_tasks/docking_run_manifest.csv
AImd/data/docking_out/file_*/*.out
AImd/data/docking_out/file_*/*_out.pdbqt
AImd/data/docking_out/docking_result_manifest.csv
AImd/data/docking_out/docking_report.json
```

`docking_result_manifest.csv` 重要字段：

```text
job_id
ligand_id
protein_id
cluster_id
batch_id
conformer_id
pocket_id
receptor_pdbqt_path
ligand_pdbqt_path
out_pose_path
log_path
best_affinity
affinities
center_x, center_y, center_z
size_x, size_y, size_z
status
```

---

## 10. Step 6：ClusterScore 提取结合能并排序重要 cluster

ClusterScore 会读取 `docking_result_manifest.csv` 或直接扫描 `.out` 文件，提取第一个构象的结合能，也就是 Vina 输出中第一行 affinity。

### 10.1 输入文件

优先读取：

```text
AImd/data/docking_out/docking_result_manifest.csv
```

如果这个文件不存在，则扫描：

```text
AImd/data/docking_out/file_*/*.out
```

### 10.2 配置文件

```text
AImd/configs/Scoring/cluster_score.yaml
```

核心参数：

```yaml
clusterscore:
  strong_threshold: -7.0
  weak_threshold: -3.0
  alpha: 0.2
  beta: 0.2
  formula: "total_strong_count + alpha * total_weak_count + n_proteins * (avg_strong_ratio + beta * avg_weak_ratio)"
  top_n: 10
```

可以改公式，例如：

```yaml
formula: "2 * total_strong_count + 0.5 * total_weak_count + n_proteins * avg_strong_ratio"
```

当前可用变量包括：

```text
alpha
beta
n_proteins
avg_strong_ratio
avg_weak_ratio
avg_none_ratio
total_strong_count
total_weak_count
total_none_count
total_binding_count
best_affinity_min
mean_affinity
```

### 10.3 运行命令

```bash
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml
```

### 10.4 主要输出

```text
AImd/data/scoring/ClusterScore/
├── best_affinity_long.csv
├── protein_ligand_matrix.csv
├── protein_binding_statistics.csv
├── cluster_binding_statistics.csv
├── top10_clusters.csv
├── clusterscore_results.xlsx
├── clusterscore_report.json
└── run_config_snapshot.json
```

其中最重要的是：

```text
clusterscore_results.xlsx
top10_clusters.csv
cluster_binding_statistics.csv
```

`clusterscore_results.xlsx` 包含：

```text
best_affinity_long
protein_ligand_matrix
protein_binding_stats
cluster_score_all
top10_clusters
```

---

## 11. Step 7：根据 Top cluster 进行精细复对接

RefinementHub 会读取 ClusterScore 的 top cluster，筛选对应蛋白，然后自动生成 refined DockingHub 配置，并再次调用 DockingHub。

### 11.1 输入文件

```text
AImd/data/scoring/ClusterScore/top10_clusters.csv
AImd/data/protein/protein_manifest.csv
AImd/data/pocket/pocket_manifest.csv
AImd/data/ligand/ligand_manifest.csv
```

### 11.2 cofactor 模板放置位置

如果精细 docking 需要 cofactor，请按 batch 放入：

```text
AImd/data/cofactor/
├── file_1/
│   ├── cyp450_heme_template_1.pdb
│   ├── fe2og_template_1.pdb
│   └── ...
├── file_2/
│   └── ...
└── ...
```

其中 `file_1`、`file_2` 要与 `protein_manifest.csv` 中的 `batch_id` 对应。

### 11.3 cofactor 模板选择逻辑

在配置中：

```yaml
cofactor:
  enabled: true
  use_foldseek: true
```

如果 `use_foldseek: true`，DockingHub 会：

```text
对当前 batch 的 cofactor 模板建 Foldseek 库
  ↓
用当前蛋白结构检索最相似模板
  ↓
用 PyMOL 对齐模板到当前蛋白
  ↓
提取模板中的 cofactor residue
  ↓
生成 receptor_with_cofactor.pdb
```

如果：

```yaml
use_foldseek: false
```

则默认使用对应 `file_*` 文件夹中的第一个模板。

### 11.4 ensemble 配置

精细 docking 默认启用 ensemble：

```yaml
ensemble:
  enabled: true
  command_template: ""
  fallback_to_input: true
  max_conformers_per_protein: 6
```

你需要根据本地 AlphaFlow 命令填写：

```yaml
command_template: "alphaflow predict --input {input_pdb} --out_dir {output_dir} --threads {threads}"
```

可用占位符：

```text
{input_pdb}
{output_dir}
{threads}
{protein_id}
{batch_id}
```

如果 `command_template` 为空，程序会使用原始结构作为 `conf_0`，不会真正产生多构象。

### 11.5 配置文件

```text
AImd/configs/Refinement/refine_from_clusterscore.yaml
```

关键参数：

```yaml
selection:
  top_n_clusters: 10
  protein_mode: all
  max_proteins_per_cluster: null

docking_overrides:
  ensemble:
    enabled: true
    align_to_reference:
      enabled: true
      fallback_to_unaligned: false

  cofactor:
    enabled: true
    use_foldseek: true

  docking:
    vina:
      exhaustiveness: 16
      cpu_per_job: 4
      timeout: 7200
```

### 11.6 运行命令

```bash
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml
```

### 11.7 主要输出

```text
AImd/data/refinement/
├── selected_clusters.csv
├── selected_protein_manifest.csv
├── refined_docking.generated.yaml
└── refinement_report.json

AImd/data/refined/
├── ensemble/conformer_manifest.csv
├── cofactor_mapped/cofactor_manifest.csv
├── receptor/receptor_manifest.csv
├── docking_configs/file_*/*.txt
├── docking_tasks/docking_task_manifest.csv
└── docking_out/docking_result_manifest.csv
```

后续 MetaBoClipBridge 会读取：

```text
AImd/data/refined/docking_out/docking_result_manifest.csv
```

---

## 12. Step 8：MetaBoClip 构象门控和家族特征打分

MetaBoClipBridge 会把 refined docking 结果整理成 MetaBoClip 的输入格式，然后运行家族特异性门控和打分。

### 12.1 输入文件

```text
AImd/data/refined/docking_out/docking_result_manifest.csv
```

其中必须包含真实存在的：

```text
receptor_pdbqt_path
out_pose_path
log_path
best_affinity
```

### 12.2 家族规则热插拔

MetaBoClip 的家族规则位于：

```text
AImd/MetaBoClip/families/
├── cyp450/
├── fe2og/
├── ugt/
├── act/
└── ach/
```

实例配置位于：

```text
AImd/MetaBoClip/configs/
├── CYP450/instance.yaml
├── Fe2OG/instance.yaml
├── UGT/instance.yaml
├── ACT/instance.yaml
└── ACH/instance.yaml
```

固定框架是：

```text
读取 receptor 和 docking pose
  ↓
识别 ligand reactive group
  ↓
识别 catalytic/cofactor atoms
  ↓
计算距离、角度、轴向、clash 等特征
  ↓
门控过滤
  ↓
pose score
  ↓
conformation score
  ↓
protein score
  ↓
final ranking
```

需要热插拔的是家族特征打分规则。新增家族时，优先新增或替换：

```text
AImd/MetaBoClip/families/<new_family>/
AImd/MetaBoClip/configs/<NEW_FAMILY>/instance.yaml
```

然后在：

```text
AImd/configs/MetaBoClip/metaboclip_bridge.yaml
```

中加入该家族。

### 12.3 配置文件

```text
AImd/configs/MetaBoClip/metaboclip_bridge.yaml
```

默认对所有内置家族都运行：

```yaml
family_assignment:
  mode: run_all

  all_families:
    - cyp450
    - fe2og
    - ugt
    - act
    - ach
```

如果你已经知道 cluster 对应的家族，可以改为：

```yaml
family_assignment:
  mode: cluster_family_map
  cluster_family_map_csv: data/input/cluster_family_map.csv
```

然后创建：

```text
AImd/data/input/cluster_family_map.csv
```

格式：

```csv
cluster_id,family
C000001,cyp450
C000002,fe2og
C000003,ugt
```

也可以一个 cluster 对应多个 family：

```csv
cluster_id,families
C000001,cyp450,fe2og
```

建议使用分号或逗号时保持一列，例如：

```csv
cluster_id,families
C000001,cyp450;fe2og
```

### 12.4 运行命令

```bash
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml
```

### 12.5 主要输出

```text
AImd/data/metaboclip/
├── staging/
│   ├── cyp450/<ligand_id>/
│   ├── fe2og/<ligand_id>/
│   └── ...
└── results/
    ├── metaboclip_run_manifest.csv
    ├── metaboclip_gating_all.csv
    ├── metaboclip_pose_scores_all.csv
    ├── metaboclip_conformation_scores_all.csv
    ├── metaboclip_protein_scores_all.csv
    ├── metaboclip_final_ranking.csv
    └── metaboclip_report.json
```

最终候选排序：

```text
AImd/data/metaboclip/results/metaboclip_final_ranking.csv
```

---

## 13. 一键运行完整流程

如果所有输入都准备好了，可以运行：

```bash
python run_full_iterative_metaboclip.py \
  --config configs/workflows/full_iterative_metaboclip.yaml
```

完整流程配置：

```text
AImd/configs/workflows/full_iterative_metaboclip.yaml
```

可控制每一步是否运行：

```yaml
workflow:
  run_rgpc: true
  run_tapocket: true
  run_broad_docking: true
  run_clusterscore: true
  run_refinement: true
  run_metaboclip: true
```

如果你已经完成前面步骤，只想从某一步继续，可以把前面的步骤改成 `false`。

例如只运行 refined docking 和 MetaBoClip：

```yaml
workflow:
  run_rgpc: false
  run_tapocket: false
  run_broad_docking: false
  run_clusterscore: false
  run_refinement: true
  run_metaboclip: true
```

---

## 14. 常用单步运行命令汇总

```bash
# 0. 检查工程结构
python validate_aimd_layout.py --root .

# 1. 蛋白结构聚类
python run_rgpc.py --config configs/RGPC/rgpc.yaml

# 2. 口袋预测
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml

# 3. 粗略分子对接
python run_docking.py --config configs/Docking/docking.yaml

# 4. ClusterScore
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml

# 5. Top cluster 精细复对接
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml

# 6. MetaBoClip 门控和家族打分
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml

# 7. 完整流程
python run_full_iterative_metaboclip.py --config configs/workflows/full_iterative_metaboclip.yaml
```

---

## 15. 最小 smoke test 建议

正式跑大规模数据前，建议准备：

```text
2-3 个蛋白结构
2-3 个 ligand PDBQT
1 个 cofactor 模板
1 个家族规则
```

先测试：

```bash
python run_rgpc.py --config configs/RGPC/rgpc.yaml
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml
python run_docking.py --config configs/Docking/docking.yaml
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml
```

确认以下文件都正常生成后，再扩大规模：

```text
data/protein/protein_manifest.csv
data/pocket/pocket_manifest.csv
data/docking_out/docking_result_manifest.csv
data/scoring/ClusterScore/clusterscore_results.xlsx
data/refined/docking_out/docking_result_manifest.csv
data/metaboclip/results/metaboclip_final_ranking.csv
```

---

## 16. 常见问题

### 16.1 运行命令后找不到数据

优先确认是否在 `AImd/` 目录下运行：

```bash
pwd
```

默认配置使用：

```yaml
aimd_root: .
```

所以工作目录必须是 AImd 根目录。

### 16.2 Foldseek 或 HipMCL 找不到

检查：

```bash
which foldseek
which hipmcl
```

如果路径不在 PATH 中，可以在配置文件里写绝对路径：

```yaml
foldseek:
  executable: /path/to/foldseek

hipmcl:
  executable: /path/to/hipmcl
```

### 16.3 Vina 没有输出 pose

检查：

```text
data/docking_out/file_*/*.out
data/docking_configs/file_*/*.txt
```

常见原因：

```text
receptor_pdbqt_path 不存在
ligand_pdbqt_path 不存在
pocket box 参数错误
vina 不在 PATH
prepare_receptor 失败
```

### 16.4 ensemble docking 后 pocket 位置不对

正确逻辑是：

```text
先用原始蛋白预测 pocket
再把 ensemble conformer 对齐回原始蛋白坐标系
然后复用原始 pocket box
```

请确认：

```yaml
ensemble:
  align_to_reference:
    enabled: true
    fallback_to_unaligned: false
```

并检查：

```text
data/refined/ensemble/conformer_manifest.csv
```

其中：

```text
alignment_status
coordinate_frame
```

应显示已经对齐到原始结构坐标系。

### 16.5 cofactor 没有映射成功

检查：

```text
data/cofactor/file_*/                 # 是否有模板
data/refined/cofactor_mapped/         # 是否有 mapped cofactor 输出
data/refined/cofactor_mapped/cofactor_manifest.csv
```

如果想暂时不因为 cofactor 映射失败而中断：

```yaml
cofactor:
  continue_without_cofactor_on_error: true
```

如果要严格要求 cofactor，改成：

```yaml
cofactor:
  continue_without_cofactor_on_error: false
```

### 16.6 ClusterScore 没有结果

确认：

```text
data/docking_out/docking_result_manifest.csv
```

中有：

```text
best_affinity
ligand_id
protein_id
cluster_id
```

如果 `best_affinity` 为空，说明 Vina `.out` 解析不到 affinity，需检查 `.out` 文件格式。

---

## 17. 后续开发建议

### 17.1 模块内部可以改，接口不要随意改

每个模块都可以内部修改：

```text
RGPC 内部换聚类算法
TApocket 内部换 pocket 预测器
DockingHub 内部换 docking engine
ScoringHub 内部换 ClusterScore 公式
MetaBoClip 内部换家族规则
```

但建议保持以下接口文件不变：

```text
data/protein/protein_manifest.csv
data/pocket/pocket_manifest.csv
data/docking_out/docking_result_manifest.csv
data/scoring/ClusterScore/top10_clusters.csv
data/refinement/selected_protein_manifest.csv
data/refined/docking_out/docking_result_manifest.csv
data/metaboclip/results/metaboclip_final_ranking.csv
```

### 17.2 工具热插拔预留

预留配置在：

```text
AImd/registry/tool_registry.yaml
AImd/registry/base.py
```

当前模块已经支持配置级替换，后续可以把工具接入统一 registry。

例如 docking engine 未来可以扩展：

```text
vina
gnina
smina
unidock
diffdock
```

pocket predictor 可以扩展：

```text
tapocket
p2rank
fpocket
pocketminer
deeppocket
```

clusterer 可以扩展：

```text
foldseek + hipmcl
foldseek easy-cluster
MMseqs2 + MCL
TM-align + Leiden
```

---

---

## 18. Third-party tool management

AImd v4 adds a centralized third-party tool layer:

```text
AImd/third_party/
├── tools.yaml
├── check_tools.py
├── setup_links.py
├── bin/
├── external/
├── resources/
└── wrappers/
```

Before running the workflow, check external tools:

```bash
cd AImd
python third_party/check_tools.py --config third_party/tools.yaml --root .
```

If a tool is installed outside `PATH`, create a stable local symlink:

```bash
python third_party/setup_links.py hipmcl /absolute/path/to/hipmcl --root .
```

Then set in `third_party/tools.yaml`:

```yaml
tools:
  hipmcl:
    executable: third_party/bin/hipmcl
```

The module configuration can use:

```yaml
executable: auto
```

When `auto` is used, AImd resolves the executable from `third_party/tools.yaml`, then from `third_party/bin`, then from system `PATH`.

## 19. RGPC Foldseek note

For RGPC, Foldseek `search` now enables alignment backtrace by default:

```yaml
foldseek:
  search:
    require_backtrace: true
```

This automatically adds `-a` to `foldseek search`. This is required by many Foldseek builds when `convertalis` exports:

```text
query,target,qtmscore,ttmscore
```

If an old `result_db` was generated without `-a`, remove old Foldseek outputs before rerunning:

```bash
rm -rf data/cluster/RGPC/foldseek_search
rm -rf data/cluster/RGPC/tmp
rm -rf data/cluster/RGPC/graph
rm -rf data/cluster/RGPC/hipmcl
python run_rgpc.py --config configs/RGPC/rgpc.yaml
```

If no non-self structural edges pass the thresholds, RGPC now skips HipMCL and creates singleton clusters automatically. This allows downstream modules to continue safely.
