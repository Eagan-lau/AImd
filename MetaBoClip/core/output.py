import pandas as pd
from pathlib import Path

def write_rows_csv(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)

def write_score_tables(score_tables, out_dir, prefix):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in score_tables.items():
        df.to_csv(out_dir / f"{prefix}.{name}.csv", index=False)
