# AlphaFlow Integration

AlphaFlow is an optional third-party conformational ensemble generator used by `DockingHub.ensemble`.
AImd does not vendor the AlphaFlow source tree or model weights. Install AlphaFlow outside the AImd
repository and point AImd to that installation through environment variables or a local config file.

Official source:

```text
https://github.com/bjing2016/alphaflow
```

## Role In AImd

AlphaFlow consumes protein sequences, MSA files, optional template PDB files, and model weights. AImd
stages these inputs from the upstream protein manifest and writes sampled conformers under:

```text
data/data_output/ensemble/
```

The stable AImd conformer manifest remains:

```text
data/data_output/ensemble/conformer_manifest.csv
```

`TApocketBridge` predicts pockets on the original protein structure. Therefore, sampled AlphaFlow
conformers are aligned back to the original protein coordinate frame before docking so that the
original pocket boxes remain valid.

## Installation

Create a dedicated AlphaFlow environment with Python 3.9:

```bash
conda create -n alphaflow python=3.9
conda activate alphaflow
```

Clone AlphaFlow outside the AImd repository:

```bash
git clone https://github.com/bjing2016/alphaflow.git /path/to/alphaflow
cd /path/to/alphaflow
```

Install dependencies following the official AlphaFlow README. The upstream project lists these core
packages:

```bash
pip install numpy==1.21.2 pandas==1.5.3
pip install torch==1.12.1+cu113 -f https://download.pytorch.org/whl/torch_stable.html
pip install biopython==1.79 dm-tree==0.1.6 modelcif==0.7 ml-collections==0.1.0 scipy==1.7.1 absl-py einops
pip install pytorch_lightning==2.0.4 fair-esm mdtraj wandb
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@103d037'
```

Download the AlphaFlow model weights into the AlphaFlow project directory. For AImd ensemble docking,
the recommended fast template-aware model is:

```text
alphaflow_12l_md_templates_base_202406.pt
```

## AImd Configuration

Export these variables before running AImd:

```bash
export ALPHAFLOW_ROOT=/path/to/alphaflow
export ALPHAFLOW_PYTHON=/path/to/miniconda/envs/alphaflow/bin/python
```

Enable the AlphaFlow adapter in `configs/Docking/docking.yaml` or a local override config:

```yaml
ensemble:
  enabled: true
  engine: alphaflow
  max_conformers_per_protein: 5
  alphaflow:
    project_dir: "${ALPHAFLOW_ROOT}"
    python_executable: "${ALPHAFLOW_PYTHON}"
    weights: alphaflow_12l_md_templates_base_202406.pt
    mode: alphafold
    samples: 5
    steps: 10
    run_msa: false
    run_prediction: true
    require_msa: true
    template_mode: copy_input
```

`run_msa: false` expects precomputed MSA files. AlphaFlow expects MSA files at:

```text
{msa_dir}/{name}/a3m/{name}.a3m
```

Set `run_msa: true` only when network access to the ColabFold MMseqs2 server is available or when
the AlphaFlow script is configured to use a local MSA workflow.

## GPU Runtime

The upstream `predict2.py` script moves the model and batches to CUDA. A compatible NVIDIA driver,
CUDA runtime, and enough GPU memory are required for real sampling. If memory is fragmented, try one
of these shell settings before running AImd:

```bash
export PYTORCH_CUDA_ALLOC_CONF=garbage_collection_threshold:0.6,max_split_size_mb:128
```

or:

```bash
export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:128
```

For multi-GPU machines:

```bash
export CUDA_VISIBLE_DEVICES=0
```

The same values can also be set in `ensemble.alphaflow.cuda_visible_devices` and
`ensemble.alphaflow.pytorch_cuda_alloc_conf`.

## Smoke Tests

Validate AImd staging without running AlphaFlow prediction:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=. python -m pytest tests/test_alphaflow_ensemble_adapter_smoke.py -q
```

Run AlphaFlow for a small local sample after MSA files and GPU support are available:

```bash
conda deactivate
conda activate aimd
python run_docking.py --config configs/Docking/docking.yaml
```

For a production workflow, use:

```bash
conda deactivate
conda activate aimd
python run_full_iterative_metaboclip.py --config configs/workflows/full_iterative_metaboclip.yaml
```

## Troubleshooting

If AlphaFlow fails with a CUDA error, confirm:

```bash
nvidia-smi
conda run -n alphaflow python -c "import torch; print(torch.cuda.is_available())"
```

If MSA files are missing, either precompute them under the configured `msa_dir` or set
`ensemble.alphaflow.run_msa: true` in an environment where the MSA query can reach its configured
server.
