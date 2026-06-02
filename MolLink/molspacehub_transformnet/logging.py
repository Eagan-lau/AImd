from __future__ import annotations

def set_rdkit_log_state(quiet: bool = True) -> None:
    """Enable/disable RDKit C++ warning logs.

    The pipeline still records parse/sanitize/key failures in QC tables. This only
    prevents console logs from hiding the actionable QC summary.
    """
    try:
        from rdkit import RDLogger
        if quiet:
            RDLogger.DisableLog("rdApp.*")
        else:
            RDLogger.EnableLog("rdApp.*")
    except Exception:
        return
