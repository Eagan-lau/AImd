from __future__ import annotations

from pathlib import Path

SCRIPT_LINES = [
    'metaboclip-build-role-table = "metaboclip_ligand_roles.build_role_table:main"',
    'metaboclip-build-role-tables-batch = "metaboclip_ligand_roles.build_role_tables_batch:main"',
    'metaboclip-recover-atom-map = "metaboclip_ligand_roles.recover_atom_map:main"',
    'metaboclip-extract-role-coords = "metaboclip_ligand_roles.extract_role_pose_coords:main"',
    'metaboclip-list-atom-labels = "metaboclip_ligand_roles.list_atom_labels:main"',
]

PACKAGE_FIND_BLOCK = [
    '[tool.setuptools.packages.find]',
    'where = ["."]',
    'include = ["metaboclip*", "metaboclip_ligand_roles*"]',
    'exclude = ["results*", "ligand_roles*", "rules*", "examples*", "docs*", "tests*", "paper_locked_originals*", "data_input*", "data_output*"]',
]


def ensure_project_scripts(text: str) -> str:
    lines = text.splitlines()
    try:
        start = lines.index('[project.scripts]')
    except ValueError:
        lines.append('')
        lines.append('[project.scripts]')
        lines.extend(SCRIPT_LINES)
        return '\n'.join(lines) + '\n'

    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].startswith('[') and lines[i].endswith(']'):
            end = i
            break

    existing_keys = set()
    for line in lines[start + 1:end]:
        if '=' in line and not line.lstrip().startswith('#'):
            existing_keys.add(line.split('=', 1)[0].strip())

    insertions = []
    for script_line in SCRIPT_LINES:
        key = script_line.split('=', 1)[0].strip()
        if key not in existing_keys:
            insertions.append(script_line)

    if insertions:
        lines[end:end] = insertions
    return '\n'.join(lines) + '\n'


def ensure_package_find(text: str) -> str:
    lines = text.splitlines()
    starts = [i for i, line in enumerate(lines) if line == '[tool.setuptools.packages.find]']
    if not starts:
        lines.append('')
        lines.extend(PACKAGE_FIND_BLOCK)
        return '\n'.join(lines) + '\n'

    start = starts[0]
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].startswith('[') and lines[i].endswith(']'):
            end = i
            break

    lines[start:end] = PACKAGE_FIND_BLOCK
    return '\n'.join(lines) + '\n'


def main() -> None:
    path = Path('pyproject.toml')
    if not path.exists():
        raise SystemExit('pyproject.toml not found. Run this script from the project root.')
    text = path.read_text(encoding='utf-8')
    text = ensure_project_scripts(text)
    text = ensure_package_find(text)
    path.write_text(text, encoding='utf-8')
    print('Updated pyproject.toml with LigandRoleMap command entries.')


if __name__ == '__main__':
    main()
