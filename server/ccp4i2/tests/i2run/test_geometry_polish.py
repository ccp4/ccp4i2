"""End-to-end geometry_polish run through i2run.

Needs Coot's headless API and MolProbity, and a refined model with its maps,
so it runs only when GEOMETRY_POLISH_E2E_DIR points at a directory holding
``model.pdb``, ``fphiout.mtz`` and (optionally) ``diffphiout.mtz``. The
Ahir3252 job 70 files (six NCS copies, 14 rotamer outliers, 6 Ramachandran
outliers) are the reference set.
"""
import csv
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from .utils import i2run

E2E_DIR = os.environ.get("GEOMETRY_POLISH_E2E_DIR")

pytestmark = pytest.mark.skipif(
    not E2E_DIR or not (Path(E2E_DIR) / "model.pdb").exists(),
    reason="set GEOMETRY_POLISH_E2E_DIR to a directory with model.pdb and fphiout.mtz")


def test_geometry_polish_end_to_end():
    data = Path(E2E_DIR)
    args = ["geometry_polish",
            "--XYZIN", str(data / "model.pdb"),
            "--FPHIOUT", str(data / "fphiout.mtz")]
    if (data / "diffphiout.mtz").exists():
        args += ["--DIFFPHIOUT", str(data / "diffphiout.mtz")]
    with i2run(args) as job:
        assert (job / "XYZOUT.pdb").exists()
        changes = job / "geometry_polish_changes.csv"
        assert changes.exists()
        with open(changes) as handle:
            rows = list(csv.DictReader(handle))
        root = ET.parse(job / "program.xml").getroot()
        before = root.find(".//summary/before")
        after = root.find(".//summary/after")
        assert before is not None and after is not None
        stages = {s.get("name"): s for s in root.findall(".//stages/stage")}
        assert set(stages) == {"rotamers", "ramachandran", "clashes"}
        # the polish may leave density-driven outliers alone, but must not add any
        assert int(after.get("rot_outliers")) <= int(before.get("rot_outliers"))
        assert int(after.get("rama_outliers")) <= int(before.get("rama_outliers"))
        accepted = root.findall(".//changes/change")
        assert len(rows) >= len(accepted)
        print("before", dict(before.attrib))
        print("after", dict(after.attrib))
        print("stages", {k: dict(v.attrib) for k, v in stages.items()})
        print("review items", len(root.findall(".//review/item")))
