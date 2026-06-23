"""
Parity test: gemmi/numpy-native CMtzData.matthewsCoeff() vs the CCP4
matthews_coef binary.

matthewsCoeff (ccp4i2.core.CCP4XtalData) reproduces matthews_coef without the
binary so the Matthews calculation runs on a CCP4-free server. matthews_coef is
the reference. For a demo MTZ's real cell/space group, across several molecular
weights and polymer modes, this asserts close agreement on Vm, percent solvent
and the Matthews probability. It runs only when the matthews_coef binary is on
PATH (i.e. under CCP4); on a slim interpreter it is skipped.
"""
import math
import os
import re
import shutil
import subprocess
import tempfile

import pytest
import gemmi

from ccp4i2 import I2_TOP
from ccp4i2.core.CCP4XtalData import CMtzDataFile

MATTHEWS = shutil.which("matthews_coef")
pytestmark = pytest.mark.skipif(MATTHEWS is None, reason="matthews_coef binary not on PATH")

MTZ = I2_TOP / "demo_data" / "gamma" / "freeR.mtz"

_ROW = re.compile(
    r'nmol_in_asu="\s*(\d+)"\s+matth_coef="\s*([0-9.E+-]+)"\s+'
    r'percent_solvent="\s*([0-9.E+-]+)"\s+prob_matth="\s*([0-9.E+-]+)"'
)


def _binary(mw, cell, sg_number, mode):
    a, b, c, al, be, ga = cell
    x = tempfile.mktemp()
    inp = (f"MOLWEIGHT {mw}\nCELL {a} {b} {c} {al} {be} {ga}\n"
           f"SYMM {sg_number}\nXMLO\nAUTO\nMODE {mode}\n")
    subprocess.run([MATTHEWS, "XMLFILE", x], input=inp, text=True,
                   capture_output=True, timeout=30)
    text = open(x).read()
    os.unlink(x)
    return [
        {"nmol_in_asu": int(m.group(1)), "matth_coef": float(m.group(2)),
         "percent_solvent": float(m.group(3)), "prob_matth": float(m.group(4))}
        for m in _ROW.finditer(text)
    ]


@pytest.fixture(scope="module")
def content():
    # Keep the CMtzDataFile alive (fileContent.cell is lost if it is GC'd).
    f = CMtzDataFile()
    f.setFullPath(str(MTZ))
    f.loadFile()
    yield f.fileContent


@pytest.mark.parametrize("mw", [4000.0, 8000.0, 16000.0, 40000.0])
@pytest.mark.parametrize("mode", ["P", "D", "C"])
def test_matthews_parity(content, mw, mode):
    mine = content.matthewsCoeff(molWt=mw, polymerMode=mode)

    cell_obj = content.cell
    angles = []
    for p in ("alpha", "beta", "gamma"):
        ang = float(getattr(cell_obj, p))
        angles.append(ang * 180.0 / math.pi if ang < 3.0 else ang)
    cell = (float(cell_obj.a), float(cell_obj.b), float(cell_obj.c), *angles)
    sg_number = gemmi.SpaceGroup(str(content.spaceGroup)).number

    ref = _binary(mw, cell, sg_number, mode)

    assert len(mine["results"]) == len(ref), "candidate count differs from binary"
    for r, b in zip(mine["results"], ref):
        assert r["nmol_in_asu"] == b["nmol_in_asu"]
        assert r["matth_coef"] == pytest.approx(b["matth_coef"], rel=1e-4)
        assert r["percent_solvent"] == pytest.approx(b["percent_solvent"], abs=2e-2)
        assert r["prob_matth"] == pytest.approx(b["prob_matth"], abs=5e-3)
