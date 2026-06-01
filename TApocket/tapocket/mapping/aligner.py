from __future__ import annotations

from pathlib import Path
from typing import Any


_PYMOL_LEGACY_LAUNCHED = False


class PyMOLAligner:
    """Map template-side objects into the original query coordinate system.

    Convention for TApocket template-only v1:
      - query is the fixed target object.
      - reference/template is the mobile object.
      - pockets/active-site objects are in the reference/template coordinate system.
      - the same transformation used for reference -> query is applied to each object.

    Runtime support:
      - Preferred: pymol2, usually provided by conda-forge pymol-open-source.
      - Fallback: legacy pymol Python module in headless mode.
    """

    def __init__(self, save_sessions: bool = False):
        self.save_sessions = save_sessions

    def _run_with_pymol2(
        self,
        query_pdb: Path,
        reference_pdb: Path,
        object_paths: list[Path],
        output_dir: Path,
        output_prefix: str,
        session_path: Path | None,
    ) -> dict[str, Any]:
        import pymol2  # type: ignore

        mapped_outputs: list[str] = []
        with pymol2.PyMOL() as pymol:
            cmd = pymol.cmd
            cmd.reinitialize()
            cmd.load(str(query_pdb), "query")
            cmd.load(str(reference_pdb), "reference")

            # Move reference into query coordinate system.
            alignment_result = cmd.super("reference", "query")
            matrix = cmd.get_object_matrix("reference")

            for index, object_path in enumerate(object_paths, start=1):
                obj_name = f"mobile_object_{index}"
                cmd.load(str(object_path), obj_name)
                cmd.transform_object(obj_name, matrix)
                out_name = f"{output_prefix}_{object_path.stem}_mapped.pdb"
                out_path = output_dir / out_name
                cmd.save(str(out_path), obj_name)
                mapped_outputs.append(str(out_path))

            if self.save_sessions and session_path:
                session_path.parent.mkdir(parents=True, exist_ok=True)
                cmd.save(str(session_path))

        return {
            "query_pdb": str(query_pdb),
            "reference_pdb": str(reference_pdb),
            "mapped_outputs": mapped_outputs,
            "alignment_result": alignment_result,
            "matrix": matrix,
            "backend": "pymol2",
        }

    def _run_with_legacy_pymol(
        self,
        query_pdb: Path,
        reference_pdb: Path,
        object_paths: list[Path],
        output_dir: Path,
        output_prefix: str,
        session_path: Path | None,
    ) -> dict[str, Any]:
        global _PYMOL_LEGACY_LAUNCHED

        import pymol  # type: ignore
        from pymol import cmd  # type: ignore

        if not _PYMOL_LEGACY_LAUNCHED:
            pymol.pymol_argv = ["pymol", "-qc"]
            pymol.finish_launching()
            _PYMOL_LEGACY_LAUNCHED = True

        mapped_outputs: list[str] = []
        cmd.reinitialize()
        cmd.load(str(query_pdb), "query")
        cmd.load(str(reference_pdb), "reference")

        # Move reference into query coordinate system.
        alignment_result = cmd.super("reference", "query")
        matrix = cmd.get_object_matrix("reference")

        for index, object_path in enumerate(object_paths, start=1):
            obj_name = f"mobile_object_{index}"
            cmd.load(str(object_path), obj_name)
            cmd.transform_object(obj_name, matrix)
            out_name = f"{output_prefix}_{object_path.stem}_mapped.pdb"
            out_path = output_dir / out_name
            cmd.save(str(out_path), obj_name)
            mapped_outputs.append(str(out_path))

        if self.save_sessions and session_path:
            session_path.parent.mkdir(parents=True, exist_ok=True)
            cmd.save(str(session_path))

        cmd.delete("all")
        return {
            "query_pdb": str(query_pdb),
            "reference_pdb": str(reference_pdb),
            "mapped_outputs": mapped_outputs,
            "alignment_result": alignment_result,
            "matrix": matrix,
            "backend": "pymol",
        }

    def map_objects_to_query(
        self,
        query_pdb: str | Path,
        reference_pdb: str | Path,
        object_paths: list[str | Path],
        output_dir: str | Path,
        output_prefix: str,
        session_path: str | Path | None = None,
    ) -> dict[str, Any]:
        query_pdb = Path(query_pdb).resolve()
        reference_pdb = Path(reference_pdb).resolve()
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        object_paths_resolved = [Path(p).resolve() for p in object_paths]
        session_path_resolved = Path(session_path).resolve() if session_path else None

        try:
            return self._run_with_pymol2(
                query_pdb=query_pdb,
                reference_pdb=reference_pdb,
                object_paths=object_paths_resolved,
                output_dir=output_dir,
                output_prefix=output_prefix,
                session_path=session_path_resolved,
            )
        except ImportError:
            pass

        try:
            return self._run_with_legacy_pymol(
                query_pdb=query_pdb,
                reference_pdb=reference_pdb,
                object_paths=object_paths_resolved,
                output_dir=output_dir,
                output_prefix=output_prefix,
                session_path=session_path_resolved,
            )
        except ImportError as exc:
            raise RuntimeError(
                "PyMOL mapping requires a Python-accessible PyMOL installation. Install it in the "
                "same environment, for example: conda install -c conda-forge pymol-open-source. "
                "Then verify with: python -c 'import pymol2; print(\"pymol2 ok\")' or "
                "python -c 'import pymol; print(\"pymol ok\")'."
            ) from exc


def get_aligner(name: str = "pymol", save_sessions: bool = False) -> PyMOLAligner:
    if name.lower() != "pymol":
        raise ValueError(f"Unsupported aligner: {name}. Currently only 'pymol' is implemented.")
    return PyMOLAligner(save_sessions=save_sessions)
