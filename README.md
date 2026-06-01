# AImd

AImd 是一个模块化的多工具蛋白-分子对接与构象打分工程包。

当前主流程：

```text
RGPC → TApocketBridge → DockingHub → ClusterScore → RefinementHub → refined DockingHub → MetaBoClipBridge
```

## 快速开始

```bash
cd AImd
pip install -r requirements.txt

python validate_aimd_layout.py --root .
```

完整运行：

```bash
python run_full_iterative_metaboclip.py \
  --config configs/workflows/full_iterative_metaboclip.yaml
```

单步运行：

```bash
python run_rgpc.py --config configs/RGPC/rgpc.yaml
python run_tapocket_batch.py --config configs/TApocket/tapocket_batch.yaml
python run_docking.py --config configs/Docking/docking.yaml
python run_clusterscore.py --config configs/Scoring/cluster_score.yaml
python run_refinement.py --config configs/Refinement/refine_from_clusterscore.yaml
python run_metaboclip_bridge.py --config configs/MetaBoClip/metaboclip_bridge.yaml
```

## 正式文档

详细使用手册：

```text
docs/AImd_USER_MANUAL.md
```

模块接口规范：

```text
docs/MODULE_INTERFACE_SPEC.md
```

工程检查报告：

```text
docs/ENGINEERING_CHECK_REPORT.md
```

## 关键输入位置

```text
data/protein_structure/cleaned_pdb/     # 蛋白结构库
data/ligand/ligand_manifest.csv         # ligand manifest
data/ligand/file_*/                     # ligand PDBQT
data/cofactor/file_*/                   # 可选 cofactor 模板
```

## 关键输出位置

```text
data/protein/protein_manifest.csv
data/pocket/pocket_manifest.csv
data/docking_out/docking_result_manifest.csv
data/scoring/ClusterScore/clusterscore_results.xlsx
data/refined/docking_out/docking_result_manifest.csv
data/metaboclip/results/metaboclip_final_ranking.csv
```

## v4 third-party layer

AImd now includes `third_party/` for external command and model management. Before running the workflow, execute:

```bash
python third_party/check_tools.py --config third_party/tools.yaml --root .
```

Use `third_party/setup_links.py` to create stable local command entries under `third_party/bin`.

RGPC now adds `-a` to Foldseek search by default (`require_backtrace: true`) and can skip HipMCL by creating singleton clusters when the filtered structural graph is empty.
