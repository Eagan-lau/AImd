# AImd 工程检查报告

本次检查对象：`AImd_Full_Iterative_MetaBoClip_integrated_v2.zip`

## 1. 清理内容

已清理：

```text
1. 顶层阶段性 README：
   - README_RGPC.md
   - README_AImd_RGPC_TApocket.md
   - README_AImd_RGPC_TApocket_Docking.md
   - README_AImd_RGPC_TApocket_Docking_ClusterScore.md
   - README_AImd_Full_Iterative_MetaBoClip.md

2. 子模块中非必要的 markdown 说明文件。
   当前正式文档统一放入 AImd/docs/。

3. Python 缓存：
   - __pycache__/
   - *.pyc

4. 构建过程遗留文件：
   - MetaBoClip/build_report.json
```

当前保留的正式文档：

```text
AImd/README.md
AImd/docs/AImd_USER_MANUAL.md
AImd/docs/MODULE_INTERFACE_SPEC.md
AImd/docs/ENGINEERING_CHECK_REPORT.md
```

## 2. 结构检查

已确认主要模块目录存在：

```text
RGPC
TApocket
TApocketBridge
DockingHub
ScoringHub
RefinementHub
MetaBoClip
MetaBoClipBridge
orchestrator
registry
configs
data
docs
```

已确认主要配置文件存在：

```text
configs/RGPC/rgpc.yaml
configs/TApocket/tapocket_batch.yaml
configs/Docking/docking.yaml
configs/Scoring/cluster_score.yaml
configs/Refinement/refine_from_clusterscore.yaml
configs/MetaBoClip/metaboclip_bridge.yaml
configs/workflows/full_iterative_metaboclip.yaml
```

## 3. 语法和导入检查

已用 Python 完成：

```text
1. 全部 .py 文件 py_compile 检查：通过
2. 全部 .yaml 文件 parse 检查：通过
3. AImd wrapper 模块导入检查：通过
```

导入检查模块：

```text
RGPC.main
TApocketBridge.main
DockingHub.main
ScoringHub.main
RefinementHub.main
MetaBoClipBridge.main
orchestrator.full_iterative_metaboclip
```

## 4. 接口 smoke test

已进行轻量接口测试：

```text
1. 构造 mock protein_manifest.csv
2. 构造 mock ligand_manifest.csv
3. 构造 mock pocket_manifest.csv
4. DockingHub 在 docking.run=false 下生成：
   - docking_task_manifest.csv
   - docking_result_manifest.csv

5. 在 docking_result_manifest.csv 中注入 mock best_affinity
6. ClusterScore 成功生成：
   - best_affinity_long.csv
   - protein_ligand_matrix.csv
   - protein_binding_statistics.csv
   - cluster_binding_statistics.csv
   - top10_clusters.csv
   - clusterscore_results.xlsx

7. RefinementHub dry-run 成功生成：
   - selected_protein_manifest.csv
   - selected_clusters.csv
   - refined_docking.generated.yaml

8. MetaBoClipBridge 在 execution.run_metaboclip=false 下成功完成 staging：
   - metaboclip_staging_manifest.csv
   - metaboclip_run_manifest.csv
   - metaboclip_final_ranking.csv
```

说明：上述 smoke test 不调用 Foldseek、HipMCL、TApocket、Vina、AlphaFlow 或 PyMOL 的真实外部计算，因此只验证模块接口和工程框架，不代表真实环境下外部工具已经可用。

## 5. 本次修正

### 5.1 ensemble/pocket 坐标规则

已固定规则：

```text
TApocket 只基于原始蛋白结构运行。
如果启用 ensemble docking，每个 conformer 必须先对齐回原始蛋白坐标系。
DockingHub 复用原始 pocket box。
```

配置中已经设置：

```yaml
ensemble:
  align_to_reference:
    enabled: true
    fallback_to_unaligned: false
```

### 5.2 防止失败 conformer 继续进入 docking

已修正 DockingHub：

```text
如果 ensemble conformer 生成失败、对齐失败或结构路径不存在，
该 conformer 不再进入 cofactor mapping / receptor preparation / docking。
```

涉及文件：

```text
DockingHub/cofactor.py
DockingHub/receptor.py
```

### 5.3 预留工具热插拔接口

新增：

```text
AImd/registry/base.py
AImd/registry/tool_registry.yaml
```

当前 workflow 仍然直接调用模块配置，registry 作为后续工具插件化的接口说明和开发起点。

## 6. 仍需在真实服务器环境确认

真实运行前必须确认：

```text
foldseek
hipmcl
pymol
prepare_receptor
vina
alphaflow   # 仅启用真实 ensemble command_template 时需要
```

还需要确认：

```text
TApocket 模板数据库和 M-CSA 数据库是否已放入对应目录；
cofactor/file_* 中是否有真实模板；
ligand_manifest.csv 中 ligand_path 是否指向真实 PDBQT；
MetaBoClip 家族 YAML 是否与当前 ligand/protein/cofactor 命名规则匹配。
```

## 7. 总体结论

当前工程包已经形成整体：

```text
RGPC → TApocketBridge → DockingHub → ClusterScore → RefinementHub → refined DockingHub → MetaBoClipBridge
```

模块之间通过 manifest 串联，后续修改可以主要限制在各模块内部。  
但真实大规模运行仍需在目标服务器上进行端到端测试，特别是外部工具路径、TApocket 数据库、AlphaFlow 命令、cofactor 模板和 MetaBoClip 家族规则。

---

## v4 update

The v4 package adds `AImd/third_party` for external tool management and patches RGPC:

- Foldseek executable can be resolved through `third_party/tools.yaml`.
- HipMCL executable can be resolved through `third_party/tools.yaml`.
- `foldseek search` adds `-a` by default via `require_backtrace: true`.
- `foldseek convertalis` defaults to `query,target,qtmscore,ttmscore`.
- RGPC parser accepts both 4-column and 5-column Foldseek TSV output.
- If the filtered similarity graph is empty, RGPC creates singleton clusters and skips HipMCL.
- DockingHub can resolve Vina, prepare_receptor, Foldseek, and PyMOL through the third-party tool registry.
