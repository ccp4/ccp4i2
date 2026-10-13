"""CRANK2 runs in CCP4i2 as the package ``crank2``, not with its own directory
on sys.path, so a bare ``from processes import x`` fails there. Its dual-mode
imports try the bare name and fall back to ``crank2.processes``; one in
processes/phdmmb.py had no fallback and every SHELX/CRANK2 SAD job whose
substructure model's space-group string differed from the data's (HypF in
H 3 2) died with "No module named 'processes'" after SHELXD had found the
sites (2026-10-13).
"""
import ast
from pathlib import Path

import ccp4i2

CRANK2 = Path(ccp4i2.__file__).parent / "pipelines" / "crank2" / "crank2"
OWN = {"processes", "programs", "basiclib"}


def _bare_own_imports(path):
    tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    guarded = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Try) and any(
                isinstance(h.type, ast.Name) and h.type.id in ("ImportError", "ModuleNotFoundError")
                for h in node.handlers):
            for inner in node.body:
                guarded.update(id(n) for n in ast.walk(inner))
    found = []
    for node in ast.walk(tree):
        if id(node) in guarded:
            continue
        if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module \
                and node.module.split(".")[0] in OWN:
            found.append(f"{path.relative_to(CRANK2)}:{node.lineno} from {node.module} import ...")
        if isinstance(node, ast.Import):
            found += [f"{path.relative_to(CRANK2)}:{node.lineno} import {a.name}"
                      for a in node.names if a.name.split(".")[0] in OWN]
    return found


def test_no_bare_import_of_crank2s_own_packages():
    bare = [hit for path in sorted(CRANK2.rglob("*.py")) for hit in _bare_own_imports(path)]
    assert not bare, "crank2 imports of its own packages need a relative or fallback form:\n" + "\n".join(bare)
